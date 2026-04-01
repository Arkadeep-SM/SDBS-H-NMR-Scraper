# SDBS ¹H NMR Scraper

A Python web scraper that automates the collection of **¹H NMR peak data** from the [SDBS (Spectral Database for Organic Compounds)](https://sdbs.db.aist.go.jp/) maintained by AIST (Japan). It extracts chemical shift (ppm), frequency (Hz), and intensity data for compounds filtered by molecular weight, then saves the results as a PyTorch tensor file and a JSON metadata file — ready for machine learning or cheminformatics workflows.

---

## Features

- Automatically accepts the SDBS disclaimer and fills the search form
- Filters compounds by **molecular weight range** (default: 200–600 Da) with **¹H NMR data available**
- Handles **multi-page pagination** across all result batches
- Navigates complex **nested frame/iframe layouts** on SDBS compound pages
- Extracts per-compound:
  - Chemical shifts (ppm), frequencies (Hz), and intensities
  - Compound name, solvent, InChI, InChIKey, and SMILES
- Saves output as:
  - `sdbs_hnmr_data.pt` — padded PyTorch tensors (hz, ppm, intensity)
  - `sdbs_hnmr_metadata.json` — structured compound metadata

---

## Requirements

- Python 3.8+
- Google Chrome (latest)

### Python Dependencies

```bash
pip install selenium webdriver-manager torch rdkit
```

Or install from a `requirements.txt`:

```bash
pip install -r requirements.txt
```

**`requirements.txt`:**
```
selenium
webdriver-manager
torch
rdkit
```

> **Note:** `rdkit` may need to be installed via conda on some platforms:
> ```bash
> conda install -c conda-forge rdkit
> ```

---

## Usage

Simply run the script:

```bash
python h1nmr_scrap_SDBS.py
```

The script will:
1. Open a Chrome browser and navigate to the SDBS disclaimer page
2. Accept the disclaimer automatically
3. Fill in the search form (MW: 200–600, ¹H NMR: checked, hits: 100 per page)
4. Iterate through all result pages and scrape peak data for each compound
5. Save the collected data to `sdbs_hnmr_data.pt` and `sdbs_hnmr_metadata.json`

### Headless Mode (Optional)

To run without a visible browser window, uncomment these lines near the top of the script:

```python
chrome_options.add_argument("--headless")
chrome_options.add_argument("--log-level=3")
```

---

## Output

### `sdbs_hnmr_data.pt`

A PyTorch `.pt` file containing a dictionary with padded tensors (zero-padded to the maximum number of peaks across all compounds):

```python
import torch

data = torch.load("sdbs_hnmr_data.pt")

data["hz"]        # torch.FloatTensor of shape [N, max_peaks]  — frequency in Hz
data["ppm"]       # torch.FloatTensor of shape [N, max_peaks]  — chemical shift in ppm
data["intensity"] # torch.IntTensor   of shape [N, max_peaks]  — relative intensity
data["metadata"]  # list of dicts, one per compound
```

### `sdbs_hnmr_metadata.json`

A JSON array of compound metadata objects:

```json
[
  {
    "sdbs_no": "12345",
    "compound_name": "Ethanol",
    "solvent": "CDCl3",
    "solvent_line": "0.5 g : 0.5 ml CDCl3",
    "inchi": "InChI=1S/C2H6O/...",
    "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    "smiles": "CCO",
    "num_peaks": 3
  },
  ...
]
```

---

## How It Works

```
SDBS Disclaimer Page
        │
        ▼
Accept Disclaimer
        │
        ▼
Fill Search Form (MW range, ¹H NMR filter, hits = 100)
        │
        ▼
Search Results Table (paginated)
        │
        ├── For each page:
        │       └── Collect SDBS numbers where HNMR = "Y"
        │               │
        │               ▼
        │          Click ¹H NMR link → Compound Page (frames/iframes)
        │               │
        │               ▼
        │          Click "peak data" button → Peak Data Page
        │               │
        │               ▼
        │          Parse: ppm, Hz, intensity, name, solvent, InChI
        │               │
        │               ▼
        │          Convert InChI → SMILES (RDKit)
        │
        ▼
Save: sdbs_hnmr_data.pt + sdbs_hnmr_metadata.json
```

---

## Configuration

You can modify these variables directly in the script:

| Parameter | Location in script | Default | Description |
|---|---|---|---|
| MW minimum | `mw_min_field.send_keys(...)` | `200` | Minimum molecular weight |
| MW maximum | `mw_max_field.send_keys(...)` | `600` | Maximum molecular weight |
| Hits per page | `Select(sel).select_by_visible_text(...)` | `100hit` | Results per search page |
| Wait timeout | `WebDriverWait(driver, 20)` | `20s` | Selenium element wait time |
| Headless mode | `chrome_options` | Disabled | Run browser without GUI |

---

## Limitations & Notes

- **Rate limiting / ToS:** This script automates access to SDBS. Please review the [SDBS Terms of Use](https://sdbs.db.aist.go.jp/Disclaimer.aspx) and use responsibly. Avoid aggressive scraping.
- **Website structure changes:** If SDBS updates its page layout, XPaths and column indices in the script may need updating.
- **Column index:** The script uses column index `5` (0-based) for the HNMR field. If results break, verify this matches the current table structure.
- **SMILES conversion:** Uses RDKit's `MolFromInchi` / `MolToSmiles`. Compounds with unparseable InChI strings will have an empty SMILES field.
- **Zero-padding:** Tensors are padded to the maximum peak count in the dataset. Use `metadata[i]["num_peaks"]` to recover the true number of peaks for compound `i`.

---

## Example: Loading and Using the Data

```python
import torch
import json

# Load tensors
data = torch.load("sdbs_hnmr_data.pt")
ppm       = data["ppm"]        # shape: [N, max_peaks]
hz        = data["hz"]         # shape: [N, max_peaks]
intensity = data["intensity"]  # shape: [N, max_peaks]
metadata  = data["metadata"]

# Load metadata
with open("sdbs_hnmr_metadata.json") as f:
    meta = json.load(f)

# Access a single compound
i = 0
n_peaks = metadata[i]["num_peaks"]
print(f"Compound: {metadata[i]['compound_name']}")
print(f"SMILES:   {metadata[i]['smiles']}")
print(f"Solvent:  {metadata[i]['solvent']}")
print(f"PPM shifts: {ppm[i, :n_peaks].tolist()}")
```

---

## Project Structure

```
.
├── h1nmr_scrap_SDBS.py         # Main scraper script
├── requirements.txt             # Python dependencies
├── README.md                    # This file
├── sdbs_hnmr_data.pt            # Output: PyTorch tensors (generated after run)
└── sdbs_hnmr_metadata.json      # Output: Compound metadata (generated after run)
```

---

## License

This project is provided for **academic and research use only**. Data collected from SDBS belongs to AIST and is subject to their terms of use. Please cite SDBS appropriately in any publications.



- [SDBS – Spectral Database for Organic Compounds](https://sdbs.db.aist.go.jp/), National Institute of Advanced Industrial Science and Technology (AIST), Japan
- [RDKit](https://www.rdkit.org/) — Open-source cheminformatics library
- [Selenium](https://www.selenium.dev/) — Browser automation framework
