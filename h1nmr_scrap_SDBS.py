from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
import time, re, json, torch
from rdkit import Chem
import pubchempy as pcp


# ── Chrome setup ───────────────────────────────────────────────────────────────
chrome_options = Options()
chrome_options.add_argument("--headless=new")
chrome_options.add_argument("--no-sandbox")
chrome_options.add_argument("--disable-dev-shm-usage")
chrome_options.add_argument("--disable-gpu")
chrome_options.add_argument("--disable-extensions")
chrome_options.add_argument("--log-level=3")
chrome_options.add_argument("--start-maximized")
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_options)
wait = WebDriverWait(driver, 20)
all_compound_data = []
# ──────────────────────────────────────────────────────────────────────────────


# ── Load existing SDBS numbers to skip duplicates ─────────────────────────────
existing_sdbs = set()
existing_data = {}
try:
    existing_data = torch.load("sdbs_hnmr_data.pt", weights_only=False)
    existing_sdbs = {str(m["sdbs_no"]) for m in existing_data["metadata"]}
    print(f"Loaded {len(existing_sdbs)} existing SDBS entries to skip.")
except FileNotFoundError:
    print("No existing .pt file found, starting fresh.")
except Exception as e:
    print(f"Warning: could not load existing .pt file: {e}")
# ──────────────────────────────────────────────────────────────────────────────


# ── PubChemPy fallback ────────────────────────────────────────────────────────
def pubchem_lookup(compound_name):
    if not compound_name:
        return None, None
    try:
        compounds = pcp.get_compounds(compound_name, "name")
        if not compounds:
            return None, None
        c = compounds[0]
        return c.canonical_smiles or None, c.inchikey or None
    except Exception as e:
        print(f"  PubChem lookup failed for '{compound_name}': {e}")
        return None, None
# ──────────────────────────────────────────────────────────────────────────────


def handle_disclaimer_if_present():
    try:
        btn = driver.find_element(By.ID, "BodyContentPlaceHolder_DisclaimeraAccept")
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
        time.sleep(0.5)
        try:
            btn.click()
        except Exception:
            driver.execute_script("arguments[0].click();", btn)
        print("Disclaimer re-accepted.")
        time.sleep(1.5)
        return True
    except Exception:
        return False


def inchi_to_smiles(inchi):
    if not inchi:
        return ""
    try:
        mol = Chem.MolFromInchi(inchi)
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception:
        return ""


def top_frame_count():
    driver.switch_to.default_content()
    return len(driver.find_elements(By.XPATH, "//frame | //iframe"))


def find_and_click_peak_data():
    handles_before = set(driver.window_handles)
    xpaths = [
        "//input[@value='peak data']",
        "//input[contains(translate(@value,'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz'),'peak')]",
        "//input[@type='submit']",
        "//button[contains(translate(text(),'ABCDEFGHIJKLMNOPQRSTUVWXYZ','abcdefghijklmnopqrstuvwxyz'),'peak')]",
    ]

    def try_click():
        for xp in xpaths:
            try:
                btn = driver.find_element(By.XPATH, xp)
                driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
                time.sleep(0.3)
                driver.execute_script("arguments[0].click();", btn)
                time.sleep(2.5)
                new_tabs = set(driver.window_handles) - handles_before
                if new_tabs:
                    driver.switch_to.window(new_tabs.pop())
                else:
                    driver.switch_to.default_content()
                return True
            except Exception:
                continue
        return False

    n = top_frame_count()
    for i in range(-1, n):
        driver.switch_to.default_content()
        if i >= 0:
            driver.switch_to.frame(i)
        if try_click():
            return True
        nested = len(driver.find_elements(By.XPATH, "//frame | //iframe"))
        for j in range(nested):
            try:
                driver.switch_to.frame(j)
                if try_click():
                    return True
                driver.switch_to.parent_frame()
            except Exception:
                driver.switch_to.default_content()
                if i >= 0:
                    driver.switch_to.frame(i)

    driver.switch_to.default_content()
    return False


def parse_peak_data_page():
    result = {}
    peak_re = re.compile(r"^\s*([\d.]+)\s+([\d.]+)\s+(\d+)\s*$", re.MULTILINE)

    n = top_frame_count()
    bodies = []
    for i in range(-1, n):
        try:
            driver.switch_to.default_content()
            if i >= 0:
                driver.switch_to.frame(i)
            bodies.append(driver.find_element(By.TAG_NAME, "body").text)
        except Exception:
            bodies.append("")
    driver.switch_to.default_content()

    best = max(bodies, key=lambda b: len(peak_re.findall(b)))
    lines = [l.strip() for l in best.split("\n") if l.strip()]

    compound_name_page = ""
    skip = re.compile(r"^(SDBS|DOI|Sample|Publisher|Rights|InChI|Hz|ppm|MHz|Assign|\d)")
    for i, line in enumerate(lines):
        if re.match(r"^[A-Z][a-z]?\d*\s+[A-Z]", line):
            for j in range(i + 1, min(i + 4, len(lines))):
                if not skip.match(lines[j]):
                    compound_name_page = lines[j]
                    break
            if compound_name_page:
                break
    result["compound_name_page"] = compound_name_page

    mhz_m = re.search(r"([\d.]+)\s*MHz", best, re.IGNORECASE)
    result["instrument_mhz"] = float(mhz_m.group(1)) if mhz_m else None

    sol_m = re.search(r"([\d.]+\s*g\s*[:\-]\s*[\d.]+\s*ml\s+(\S+))", best)
    if sol_m:
        result["solvent_line"] = sol_m.group(1).strip()
        result["solvent"] = sol_m.group(2).strip()
    else:
        s = re.search(
            r"\b(CDCl3|D2O|DMSO(?:-d6)?|CD3OD|C6D6|CD2Cl2|acetone-d6|methanol-d4)\b",
            best, re.IGNORECASE,
        )
        result["solvent"] = s.group(1) if s else "unknown"
        result["solvent_line"] = result["solvent"]

    hz_list, ppm_list, int_list = [], [], []
    for m in peak_re.finditer(best):
        hz_list.append(float(m.group(1)))
        ppm_list.append(float(m.group(2)))
        int_list.append(int(m.group(3)))
    result["hz"], result["ppm"], result["int"] = hz_list, ppm_list, int_list

    inchi_m = re.search(r"(InChI=1S?/[^\s\n]+)", best)
    result["inchi"] = inchi_m.group(1).strip() if inchi_m else ""

    ikey_m = re.search(r"([A-Z]{14}-[A-Z]{10}-[A-Z])", best)
    result["inchikey"] = ikey_m.group(1).strip() if ikey_m else ""

    return result


def scrape_hnmr_sdbs_numbers():
    """
    Returns list of dicts: {"sdbs_no": str, "compound_name": str}
    for every row where HNMR column == "Y".

    Column layout:
      0  SDBS No
      1  Molecular Formula
      2  Molecular Weight
      3  MS
      4  CNMR
      5  HNMR   <- Y/N
      6  IR
      7  Raman
      8  ESR
      9  Compound Name
    """
    sdbs_list = []
    try:
        wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
        time.sleep(1)
        rows = driver.find_elements(By.XPATH, "//table//tr[td]")
        for row in rows:
            cols = row.find_elements(By.TAG_NAME, "td")
            if len(cols) >= 6:
                sdbs_no = cols[0].text.strip()
                hnmr_val = cols[5].text.strip().upper()
                compound_name = cols[9].text.strip() if len(cols) > 9 else ""
                if hnmr_val == "Y" and sdbs_no.isdigit():
                    sdbs_list.append({"sdbs_no": sdbs_no, "compound_name": compound_name})
    except Exception as e:
        print(f"scrape error: {e}")
    return sdbs_list


def click_hnmr_y_for_sdbs(sdbs_no):
    try:
        wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
        rows = driver.find_elements(By.XPATH, "//table//tr[td]")
        for row in rows:
            cols = row.find_elements(By.TAG_NAME, "td")
            if len(cols) >= 6 and cols[0].text.strip() == sdbs_no:
                hnmr_col = cols[5]
                try:
                    link = hnmr_col.find_element(By.TAG_NAME, "a")
                    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", link)
                    time.sleep(0.2)
                    driver.execute_script("arguments[0].click();", link)
                except Exception:
                    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", hnmr_col)
                    time.sleep(0.2)
                    driver.execute_script("arguments[0].click();", hnmr_col)
                return True
    except Exception as e:
        print(f"click_hnmr_y({sdbs_no}): {e}")
    return False


# ── Pagination helpers (numbered pages) ───────────────────────────────────────

def get_all_page_numbers():
    """
    Reads the pagination bar and returns a sorted list of all integer page
    numbers (1-based), including the current page which is NOT a link.

    The bar looks like:  1  2  3  4  5  6
    where 1 is the current page (plain text / span, NOT an <a>) and 2-6
    are <a href="..."> links.

    Key insight: pagination <a> links point back to SearchResult.aspx,
    while compound-row <a> links always contain 'sdbsno' or 'Details' in
    their href.  Filter on href to exclude table links entirely.
    Then infer the current (non-linked) page number from JS text-node scan
    or from the gap in the sequence.
    """
    try:
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(0.4)

        page_nums = set()

        # ── Step 1: collect pagination <a> links by href pattern ───────────
        # Compound links:   href contains 'sdbsno' or 'Details'
        # Pagination links: href contains 'SearchResult' OR is a JS postback
        #                   (javascript:__doPostBack) — neither has 'sdbsno'
        all_links = driver.find_elements(By.TAG_NAME, "a")
        for a in all_links:
            txt = a.text.strip()
            if not txt.isdigit():
                continue
            num = int(txt)
            if num < 1 or num > 500:          # page numbers are always small
                continue
            href = (a.get_attribute("href") or "").lower()
            # Hard exclude: any link that goes to a compound detail page
            if "sdbsno" in href or "details" in href:
                continue
            page_nums.add(num)

        # ── Step 2: find the current page number (not a link) via JS ───────
        # Try several common patterns: <b>, <strong>, <span>, bare text node
        current_page = driver.execute_script("""
            // 1. Check <b> / <strong> / <span> with only digit content
            var tags = ['b','strong','span','font','td'];
            for (var t of tags) {
                var els = document.getElementsByTagName(t);
                for (var el of els) {
                    var txt = el.textContent.trim();
                    if (/^\\d+$/.test(txt) && parseInt(txt) >= 1
                            && parseInt(txt) <= 500
                            && el.tagName.toLowerCase() !== 'a') {
                        // Make sure it has no <a> child (i.e. it IS the current page)
                        if (!el.querySelector('a')) {
                            // Heuristic: must be near other small-digit siblings/neighbours
                            var parent = el.parentNode;
                            if (parent) {
                                var siblings = parent.querySelectorAll('a');
                                var hasPageLinks = false;
                                for (var s of siblings) {
                                    var st = s.textContent.trim();
                                    if (/^\\d+$/.test(st) && parseInt(st) <= 500) {
                                        hasPageLinks = true; break;
                                    }
                                }
                                if (hasPageLinks) return parseInt(txt);
                            }
                        }
                    }
                }
            }
            // 2. Walk text nodes near the bottom of the body
            var walker = document.createTreeWalker(
                document.body, NodeFilter.SHOW_TEXT, null, false);
            var node;
            while (node = walker.nextNode()) {
                var txt = node.textContent.trim();
                if (/^\\d+$/.test(txt)) {
                    var n = parseInt(txt);
                    if (n >= 1 && n <= 500) {
                        var p = node.parentNode;
                        if (p && p.tagName && p.tagName.toLowerCase() !== 'a') {
                            return n;
                        }
                    }
                }
            }
            return null;
        """)

        if current_page and isinstance(current_page, (int, float)):
            page_nums.add(int(current_page))

        # ── Step 3: if we found links (e.g. 2-6) but no current page, ──────
        #           the current page is 1 (we always start from page 1)
        if page_nums and 1 not in page_nums:
            # Check whether the sequence starts at 2 → current page is 1
            if min(page_nums) == 2:
                page_nums.add(1)

        if not page_nums:
            print("  [pagination] Could not detect page numbers; assuming 1 page.")
            return [1]

        return sorted(page_nums)

    except Exception as e:
        print(f"get_all_page_numbers error: {e}")
        return [1]


def get_current_page_number():
    """
    Returns the currently active page number by looking for the bold /
    non-linked digit in the pagination bar.  Falls back to 1.
    """
    try:
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(0.3)

        # Try common patterns: <b>, <strong>, <span class="...current...">
        for xpath in [
            "//b[number(normalize-space(text())) = number(normalize-space(text()))]",
            "//strong[number(normalize-space(text())) = number(normalize-space(text()))]",
            "//*[contains(@class,'current') or contains(@class,'active') or contains(@class,'selected')]"
            "[number(normalize-space(text())) = number(normalize-space(text()))]",
        ]:
            els = driver.find_elements(By.XPATH, xpath)
            for el in els:
                txt = el.text.strip()
                if txt.isdigit():
                    return int(txt)
    except Exception:
        pass
    return 1


def click_page_number(page_num):
    """
    Clicks the <a> link whose text is exactly `page_num`.
    Returns True on success, False otherwise.
    """
    try:
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(0.4)
        link = driver.find_element(
            By.XPATH,
            f"//a[normalize-space(text())='{page_num}']"
        )
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", link)
        time.sleep(0.3)
        driver.execute_script("arguments[0].click();", link)
        time.sleep(2.5)
        handle_disclaimer_if_present()
        wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
        time.sleep(1)
        return True
    except Exception as e:
        print(f"click_page_number({page_num}) error: {e}")
        return False

# ──────────────────────────────────────────────────────────────────────────────


def process_compound(sdbs_no, compound_name_from_table, results_url):
    main_handle = driver.window_handles[0]
    try:
        handle_disclaimer_if_present()

        if driver.current_url != results_url:
            driver.get(results_url)
            time.sleep(1.5)
            handle_disclaimer_if_present()
            wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
            time.sleep(1.5)

        handles_before_comp = set(driver.window_handles)
        if not click_hnmr_y_for_sdbs(sdbs_no):
            print(f"  [{sdbs_no}] HNMR Y link not found")
            return None
        time.sleep(2.5)

        new_comp = set(driver.window_handles) - handles_before_comp
        if new_comp:
            driver.switch_to.window(new_comp.pop())

        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        time.sleep(2)
        handle_disclaimer_if_present()

        if not find_and_click_peak_data():
            print(f"  [{sdbs_no}] peak data button not found")
            return None

        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        time.sleep(1.5)
        handle_disclaimer_if_present()

        data = parse_peak_data_page()
        data["sdbs_no"] = sdbs_no
        data["compound_name"] = compound_name_from_table or data.pop("compound_name_page", "")

        if not data.get("inchikey"):
            print(f"  [{sdbs_no}] InChIKey missing — querying PubChem for '{data['compound_name']}'...")
            pc_smiles, pc_inchikey = pubchem_lookup(data["compound_name"])
            if pc_inchikey:
                data["inchikey"] = pc_inchikey
                print(f"  [{sdbs_no}] PubChem InChIKey: {pc_inchikey}")
            if pc_smiles and not data.get("inchi"):
                data["smiles_pubchem"] = pc_smiles
                print(f"  [{sdbs_no}] PubChem SMILES:   {pc_smiles}")

        print(f"  [{sdbs_no}] {data.get('compound_name','?')} | "
              f"{len(data.get('ppm', []))} peaks | "
              f"{data.get('solvent','?')} | "
              f"{data.get('instrument_mhz','?')} MHz")
        return data

    except Exception as e:
        print(f"  [{sdbs_no}] error: {e}")
        return None

    finally:
        for handle in list(driver.window_handles):
            if handle != main_handle:
                driver.switch_to.window(handle)
                driver.close()
        driver.switch_to.window(main_handle)
        try:
            if driver.current_url != results_url:
                driver.get(results_url)
                time.sleep(1.5)
                handle_disclaimer_if_present()
                wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
                time.sleep(1.5)
        except Exception:
            pass


# ── Main ───────────────────────────────────────────────────────────────────────
try:
    print("Navigating to SDBS...")
    driver.get("https://sdbs.db.aist.go.jp/Disclaimer.aspx")
    button = wait.until(EC.presence_of_element_located(
        (By.ID, "BodyContentPlaceHolder_DisclaimeraAccept")
    ))
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    time.sleep(1)
    try:
        button.click()
    except Exception:
        driver.execute_script("arguments[0].click();", button)
    print("Disclaimer accepted.")
    time.sleep(2)

    # ── Fill search form ───────────────────────────────────────────────────────
    all_inputs = driver.find_elements(By.XPATH, "//input[@type='text' or @type='']")
    mw_min_field = mw_max_field = None
    for inp in all_inputs:
        val  = inp.get_attribute("value") or ""
        name = (inp.get_attribute("name") or "").lower()
        id_  = (inp.get_attribute("id") or "").lower()
        if "200" in val or "min" in name or "min" in id_:
            mw_min_field = inp
        if "900" in val or "max" in name or "max" in id_:
            mw_max_field = inp
    if not mw_min_field or not mw_max_field:
        mw_label  = driver.find_element(By.XPATH, "//*[contains(text(),'Molecular Weight')]")
        mw_inputs = mw_label.find_elements(
            By.XPATH, "following::input[@type='text'][position()<=2]"
        )
        if len(mw_inputs) >= 2:
            mw_min_field, mw_max_field = mw_inputs[0], mw_inputs[1]
    if mw_min_field:
        mw_min_field.clear(); mw_min_field.send_keys("200")
    if mw_max_field:
        mw_max_field.clear(); mw_max_field.send_keys("900")
    print("MW 200-900 set.")

    checkboxes = driver.find_elements(By.XPATH, "//input[@type='checkbox']")
    hnmr_cb = None
    for cb in checkboxes:
        combined = (
            (cb.get_attribute("id") or "")
            + (cb.get_attribute("name") or "")
            + (cb.get_attribute("value") or "")
        ).lower()
        if any(k in combined for k in ["hnmr", "1hnmr", "1h_nmr", "h1nmr", "proton"]):
            hnmr_cb = cb
            break
    if not hnmr_cb and len(checkboxes) >= 5:
        hnmr_cb = checkboxes[4]
    if hnmr_cb and not hnmr_cb.is_selected():
        driver.execute_script("arguments[0].click();", hnmr_cb)
    print("1H NMR checked.")

    for sel in driver.find_elements(By.TAG_NAME, "select"):
        opts = [o.text for o in Select(sel).options]
        if any("100" in o for o in opts):
            Select(sel).select_by_visible_text("100hit")
            break
    print("Hit = 100hit.")

    # ── Click Search ───────────────────────────────────────────────────────────
    tabs_before = set(driver.window_handles)
    url_before  = driver.current_url
    search_btn  = driver.find_element(
        By.XPATH,
        "//input[@type='submit' and @value='Search'] | "
        "//button[contains(text(),'Search')] | "
        "//input[@type='button' and @value='Search']",
    )
    driver.execute_script("arguments[0].click();", search_btn)
    print("Search clicked...")
    for _ in range(30):
        time.sleep(1)
        new_tabs = set(driver.window_handles) - tabs_before
        if new_tabs:
            driver.switch_to.window(new_tabs.pop())
            break
        if driver.current_url != url_before:
            break
    time.sleep(2)
    handle_disclaimer_if_present()
    print(f"Results page: {driver.current_url}")
    # ──────────────────────────────────────────────────────────────────────────

    # ── Pagination loop — page numbers re-detected fresh on every page ────────
    current_page_num = 1

    while True:

        # ── Detect pages from whatever is currently loaded ─────────────────
        all_pages   = get_all_page_numbers() or [1]
        total_pages = max(all_pages)
        remaining   = [p for p in all_pages if p > current_page_num]

        print(f"\n{'='*55}")
        print(f"  Page numbers detected : {all_pages}")
        print(f"  Currently on page     : {current_page_num} / {total_pages}")
        if remaining:
            print(f"  Pages still to visit  : {remaining}")
        else:
            print(f"  Pages still to visit  : none — this is the last page")
        print(f"{'='*55}")
        # ──────────────────────────────────────────────────────────────────

        current_url = driver.current_url
        print(f"\n=== Processing Page {current_page_num}/{total_pages} | URL: {current_url} ===")

        sdbs_list = scrape_hnmr_sdbs_numbers()
        print(f"  Found {len(sdbs_list)} HNMR=Y entries on page {current_page_num}")

        for i, entry in enumerate(sdbs_list):
            sdbs_no       = entry["sdbs_no"]
            compound_name = entry["compound_name"]

            if sdbs_no in existing_sdbs:
                print(f"  [{i+1}/{len(sdbs_list)}] SDBS {sdbs_no} — skipped (already exists)")
                continue

            print(f"  [{i+1}/{len(sdbs_list)}] SDBS {sdbs_no} — {compound_name}")
            data = process_compound(sdbs_no, compound_name, current_url)
            if data:
                all_compound_data.append(data)
                existing_sdbs.add(sdbs_no)

        # ── Decide where to go next ────────────────────────────────────────
        if not remaining:
            print(f"\n  No more pages. All done.")
            break

        next_page = remaining[0]   # always the immediate next page number
        print(f"\n{'─'*45}")
        print(f"  Navigating to page {next_page} of {total_pages}...")
        print(f"  Pages left after this : {remaining[1:] if len(remaining) > 1 else 'none'}")
        print(f"{'─'*45}")

        if not click_page_number(next_page):
            print(f"  Could not click page {next_page}. Stopping.")
            break

        current_page_num = next_page
        # ──────────────────────────────────────────────────────────────────

    print(f"\n{'='*55}")
    print(f"  Finished all {total_pages} pages.")
    print(f"  Total new compounds scraped: {len(all_compound_data)}")
    print(f"{'='*55}")
    # ──────────────────────────────────────────────────────────────────────────

    # ── Merge with existing data ───────────────────────────────────────────────
    if existing_data and existing_data.get("metadata"):
        print(f"Merging with {len(existing_data['metadata'])} existing entries...")
        ex_meta = existing_data["metadata"]
        ex_hz   = existing_data["hz"].tolist()
        ex_ppm  = existing_data["ppm"].tolist()
        ex_int  = existing_data["intensity"].tolist()
        merged = []
        for j, m in enumerate(ex_meta):
            n = m["num_peaks"]
            merged.append({
                "sdbs_no":        m["sdbs_no"],
                "compound_name":  m["compound_name"],
                "solvent":        m["solvent"],
                "solvent_line":   m["solvent_line"],
                "inchi":          m["inchi"],
                "inchikey":       m["inchikey"],
                "instrument_mhz": m.get("instrument_mhz"),
                "ppm":            ex_ppm[j][:n],
                "hz":             ex_hz[j][:n],
                "int":            ex_int[j][:n],
            })
        merged.extend(all_compound_data)
        all_compound_data_final = merged
    else:
        all_compound_data_final = all_compound_data
    # ──────────────────────────────────────────────────────────────────────────

    # ── Save ───────────────────────────────────────────────────────────────────
    print(f"Total compounds (merged): {len(all_compound_data_final)}")
    max_peaks = max((len(d.get("ppm", [])) for d in all_compound_data_final), default=0)
    hz_t, ppm_t, int_t, metadata = [], [], [], []
    for d in all_compound_data_final:
        n   = len(d.get("ppm", []))
        pad = max_peaks - n
        hz_t.append(d.get("hz",  []) + [0.0] * pad)
        ppm_t.append(d.get("ppm", []) + [0.0] * pad)
        int_t.append(d.get("int", []) + [0]   * pad)

        smiles = inchi_to_smiles(d.get("inchi", "")) or d.get("smiles_pubchem", "")

        metadata.append({
            "sdbs_no":        d.get("sdbs_no", ""),
            "compound_name":  d.get("compound_name", ""),
            "solvent":        d.get("solvent", ""),
            "solvent_line":   d.get("solvent_line", ""),
            "inchi":          d.get("inchi", ""),
            "inchikey":       d.get("inchikey", ""),
            "smiles":         smiles,
            "instrument_mhz": d.get("instrument_mhz"),
            "num_peaks":      n,
        })

    torch.save({
        "hz":        torch.tensor(hz_t,  dtype=torch.float32),
        "ppm":       torch.tensor(ppm_t, dtype=torch.float32),
        "intensity": torch.tensor(int_t, dtype=torch.int32),
        "metadata":  metadata,
    }, "temp_sdbs_hnmr_data.pt")

    with open("temp_sdbs_hnmr_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    print("Saved: temp_sdbs_hnmr_data.pt  |  temp_sdbs_hnmr_metadata.json")

except Exception as e:
    import traceback
    print(f"Fatal: {e}")
    traceback.print_exc()

finally:
    driver.quit()
    print("Browser closed.")
