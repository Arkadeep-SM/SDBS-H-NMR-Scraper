from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
import time, re, json, torch
from rdkit import Chem


chrome_options = Options()
#chrome_options.add_argument("--headless")
#chrome_options.add_argument("--log-level=3")  # Set the log level to 'warning'.
chrome_options.add_argument('--start-maximized')
driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_options)
wait = WebDriverWait(driver, 20)
all_compound_data = []





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
    peak_re = re.compile(r'^\s*([\d.]+)\s+([\d.]+)\s+(\d+)\s*$', re.MULTILINE)

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

    compound_name = ""
    skip = re.compile(r'^(SDBS|DOI|Sample|Publisher|Rights|InChI|Hz|ppm|MHz|Assign|\d)')
    for i, line in enumerate(lines):
        if re.match(r'^[A-Z][a-z]?\d*\s+[A-Z]', line):
            for j in range(i + 1, min(i + 4, len(lines))):
                if not skip.match(lines[j]):
                    compound_name = lines[j]
                    break
            if compound_name:
                break
    result["compound_name"] = compound_name

    sol_m = re.search(r'([\d.]+\s*g\s*[:\-]\s*[\d.]+\s*ml\s+(\S+))', best)
    if sol_m:
        result["solvent_line"] = sol_m.group(1).strip()
        result["solvent"]      = sol_m.group(2).strip()
    else:
        s = re.search(
            r'\b(CDCl3|D2O|DMSO(?:-d6)?|CD3OD|C6D6|CD2Cl2|acetone-d6|methanol-d4)\b',
            best, re.IGNORECASE
        )
        result["solvent"]      = s.group(1) if s else "unknown"
        result["solvent_line"] = result["solvent"]

    hz_list, ppm_list, int_list = [], [], []
    for m in peak_re.finditer(best):
        hz_list.append(float(m.group(1)))
        ppm_list.append(float(m.group(2)))
        int_list.append(int(m.group(3)))
    result["hz"], result["ppm"], result["int"] = hz_list, ppm_list, int_list

    inchi_m = re.search(r'(InChI=1S?/[^\s\n]+)', best)
    result["inchi"] = inchi_m.group(1).strip() if inchi_m else ""

    ikey_m = re.search(r'([A-Z]{14}-[A-Z]{10}-[A-Z])', best)
    result["inchikey"] = ikey_m.group(1).strip() if ikey_m else ""

    return result


def scrape_hnmr_sdbs_numbers():
    sdbs_list = []
    try:
        wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
        time.sleep(1)
        rows = driver.find_elements(By.XPATH, "//table//tr[td]")
        for row in rows:
            cols = row.find_elements(By.TAG_NAME, "td")
            if len(cols) >= 6:
                sdbs_no  = cols[0].text.strip()
                hnmr_val = cols[5].text.strip().upper()  # was cols[4]
                if hnmr_val == "Y" and sdbs_no.isdigit():
                    sdbs_list.append(sdbs_no)
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
                hnmr_col = cols[5]  # was cols[4]
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

def get_pagination_info():
    try:
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(0.5)
        page_links = driver.find_elements(
            By.XPATH,
            "//a[string-length(normalize-space(.)) <= 4 and "
            "translate(normalize-space(.), '0123456789', '') = '']"
        )
        visible_pages = sorted(
            set(int(a.text.strip()) for a in page_links if a.text.strip().isdigit())
        )
        nxt = driver.find_elements(
            By.XPATH,
            "//a[normalize-space(text())='»' or "
            "normalize-space(text())='>>' or "
            "normalize-space(text())='Next']"
        )
        has_next = len(nxt) > 0
        return visible_pages, has_next, (nxt[0] if has_next else None)
    except Exception:
        return [], False, None


def click_page_number(page_num):
    try:
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(0.3)
        link = driver.find_element(
            By.XPATH,
            f"//a[normalize-space(text())='{page_num}' and "
            f"translate(normalize-space(text()), '0123456789', '') = '']"
        )
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", link)
        time.sleep(0.3)
        driver.execute_script("arguments[0].click();", link)
        time.sleep(2)
        return True
    except Exception as e:
        print(f"click_page({page_num}): {e}")
        return False


def process_compound(sdbs_no, results_url):
    main_handle = driver.window_handles[0]
    try:
        if driver.current_url != results_url:
            driver.get(results_url)
            wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
            time.sleep(1.5)

        handles_before_comp = set(driver.window_handles)
        if not click_hnmr_y_for_sdbs(sdbs_no):
            print(f"[{sdbs_no}] HNMR Y link not found")
            return None
        time.sleep(2.5)

        new_comp = set(driver.window_handles) - handles_before_comp
        if new_comp:
            driver.switch_to.window(new_comp.pop())

        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        time.sleep(2)

        if not find_and_click_peak_data():
            print(f"[{sdbs_no}] peak data button not found")
            return None

        wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
        time.sleep(1.5)

        data = parse_peak_data_page()
        data["sdbs_no"] = sdbs_no
        print(f"[{sdbs_no}] {data.get('compound_name','?')} | {len(data.get('ppm',[]))} peaks | {data.get('solvent','?')}")
        return data

    except Exception as e:
        print(f"[{sdbs_no}] error: {e}")
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
                wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
                time.sleep(1.5)
        except Exception:
            pass


try:
    # Disclaimer
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

    # Fill form
    all_inputs = driver.find_elements(By.XPATH, "//input[@type='text' or @type='']")
    mw_min_field = mw_max_field = None
    for inp in all_inputs:
        val  = inp.get_attribute('value') or ''
        name = (inp.get_attribute('name') or '').lower()
        id_  = (inp.get_attribute('id') or '').lower()
        if '200' in val or 'min' in name or 'min' in id_:
            mw_min_field = inp
        if '600' in val or 'max' in name or 'max' in id_:
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
        mw_max_field.clear(); mw_max_field.send_keys("600")
    print("MW 200-600 set.")

    checkboxes = driver.find_elements(By.XPATH, "//input[@type='checkbox']")
    hnmr_cb = None
    for cb in checkboxes:
        combined = ((cb.get_attribute('id') or '') +
                    (cb.get_attribute('name') or '') +
                    (cb.get_attribute('value') or '')).lower()
        if any(k in combined for k in ['hnmr','1hnmr','1h_nmr','h1nmr','proton']):
            hnmr_cb = cb; break
    if not hnmr_cb and len(checkboxes) >= 5:
        hnmr_cb = checkboxes[4]
    if hnmr_cb and not hnmr_cb.is_selected():
        driver.execute_script("arguments[0].click();", hnmr_cb)
    print("1H NMR checked.")

    for sel in driver.find_elements(By.TAG_NAME, "select"):
        opts = [o.text for o in Select(sel).options]
        if any('100' in o for o in opts):
            Select(sel).select_by_visible_text("100hit"); break
    print("Hit = 100hit.")

    # Search
    tabs_before = set(driver.window_handles)
    url_before  = driver.current_url
    search_btn  = driver.find_element(
        By.XPATH,
        "//input[@type='submit' and @value='Search'] | "
        "//button[contains(text(),'Search')] | "
        "//input[@type='button' and @value='Search']"
    )
    driver.execute_script("arguments[0].click();", search_btn)
    print("Search clicked...")
    for _ in range(30):
        time.sleep(1)
        new_tabs = set(driver.window_handles) - tabs_before
        if new_tabs:
            driver.switch_to.window(new_tabs.pop()); break
        if driver.current_url != url_before:
            break
    time.sleep(2)
    results_url = driver.current_url
    print(f"Results: {results_url}")

    # Paginate and process
    global_page = 0
    while True:
        visible_pages, has_next_batch, next_batch_el = get_pagination_info()

        if not visible_pages:
            print("Single page.")
            sdbs_list = scrape_hnmr_sdbs_numbers()
            print(f"{len(sdbs_list)} HNMR=Y entries")
            for i, sdbs_no in enumerate(sdbs_list):
                print(f"[{i+1}/{len(sdbs_list)}] SDBS {sdbs_no}")
                data = process_compound(sdbs_no, results_url)
                if data:
                    all_compound_data.append(data)
            break

        print(f"Batch: {visible_pages} | has_next: {has_next_batch}")

        for p in visible_pages:
            global_page += 1
            print(f"Page {p} (#{global_page})")
            if p != visible_pages[0] or global_page > 1:
                if not click_page_number(p):
                    continue
                wait.until(EC.presence_of_element_located((By.XPATH, "//table//tr[td]")))
                time.sleep(1)
            results_url = driver.current_url
            sdbs_list   = scrape_hnmr_sdbs_numbers()
            print(f"{len(sdbs_list)} HNMR=Y: {sdbs_list}")
            for i, sdbs_no in enumerate(sdbs_list):
                print(f"[{i+1}/{len(sdbs_list)}] SDBS {sdbs_no}")
                data = process_compound(sdbs_no, results_url)
                if data:
                    all_compound_data.append(data)

        if has_next_batch:
            print("Next batch...")
            try:
                driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
                time.sleep(0.5)
                nxt = driver.find_element(
                    By.XPATH,
                    "//a[normalize-space(text())='»' or "
                    "normalize-space(text())='>>' or "
                    "normalize-space(text())='Next']"
                )
                driver.execute_script("arguments[0].click();", nxt)
                time.sleep(2)
                results_url = driver.current_url
            except Exception as e:
                print(f"next batch error: {e}"); break
        else:
            print("All done.")
            break

    # Save
    print(f"Total compounds: {len(all_compound_data)}")
    max_peaks = max((len(d.get("ppm", [])) for d in all_compound_data), default=0)
    hz_t, ppm_t, int_t, metadata = [], [], [], []
    for d in all_compound_data:
        n   = len(d.get("ppm", []))
        pad = max_peaks - n
        hz_t.append(d.get("hz",  []) + [0.0] * pad)
        ppm_t.append(d.get("ppm", []) + [0.0] * pad)
        int_t.append(d.get("int", []) + [0]   * pad)
        metadata.append({
            "sdbs_no"      : d.get("sdbs_no", ""),
            "compound_name": d.get("compound_name", ""),
            "solvent"      : d.get("solvent", ""),
            "solvent_line" : d.get("solvent_line", ""),
            "inchi"        : d.get("inchi", ""),
            "inchikey"     : d.get("inchikey", ""),
            "smiles"       : inchi_to_smiles(d.get("inchi", "")),  # <-- add this
            "num_peaks"    : n,
        })

    torch.save({
        "hz"       : torch.tensor(hz_t,  dtype=torch.float32),
        "ppm"      : torch.tensor(ppm_t, dtype=torch.float32),
        "intensity": torch.tensor(int_t, dtype=torch.int32),
        "metadata" : metadata,
    }, "sdbs_hnmr_data.pt")

    with open("sdbs_hnmr_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    print("Saved: sdbs_hnmr_data.pt, sdbs_hnmr_metadata.json")

except Exception as e:
    import traceback
    print(f"Fatal: {e}")
    traceback.print_exc()

finally:
    driver.quit()
    print("Browser closed.")