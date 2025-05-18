# Directory Structures

This document describes the directory layouts as they appear on PhysioNet for three MIMIC waveform datasets:

---

## MIMIC-III Waveform Database Matched Subset v1.0

PhysioNet path: `https://physionet.org/files/mimic3wdb-matched/1.0/` ([archive.physionet.org](https://archive.physionet.org/physiobank/database/mimic3wdb/matched/?utm_source=chatgpt.com))

```plaintext
mimic3wdb-matched/
└── 1.0/
    ├── p00/                # Intermediate directories p00–p09
    │   ├── p00NNNN/        # Subject directories (Subject_ID = XXNNNN)
    │   │   ├── p00NNNN-YYYY-MM-DD-hh-mm.hea  # WFDB header file
    │   │   ├── p00NNNN-YYYY-MM-DD-hh-mm.dat  # WFDB signal data
    │   │   └── p00NNNN-YYYY-MM-DD-hh-mm.n    # Numerics record (suffix 'n')
    │   └── ...
    └── p01/ … p09/
```

Notes:

* Records are organized in ten `pXX/` folders to reduce directory sizes.
* Each dated record has both `.hea` and `.dat` files, plus an optional numerics `.n` file. ([archive.physionet.org](https://archive.physionet.org/physiobank/database/mimic3wdb/matched/?utm_source=chatgpt.com))

---

## MIMIC-IV Waveform Database v0.1.0

PhysioNet path: `https://physionet.org/files/mimic4wdb/0.1.0/` ([physionet.org](https://www.physionet.org/content/mimic4wdb/0.1.0/waves/p100/p10039708/83411188/83411188_0001e.dat?utm_source=chatgpt.com))

```plaintext
mimic4wdb/
└── 0.1.0/
    └── waves/
        ├── p100/               # First-level grouping (e.g., by record prefix)
        │   ├── p10039708/      # Secondary-level directory
        │   │   ├── 83411188/   # Record folder (record number)
        │   │   │   ├── 83411188_0001e.dat  # Signal data file
        │   │   │   └── 83411188_0001e.hea  # Header file
        │   │   └── ...
        │   └── ...
        └── ...
```

Notes:

* The `waves/` directory holds WFDB-formatted signals, grouped hierarchically for performance. ([physionet.org](https://www.physionet.org/content/mimic4wdb/0.1.0/waves/p100/p10039708/83411188/83411188_0001e.dat?utm_source=chatgpt.com))

---

## MIMIC-IV ECG: Diagnostic ECG Matched Subset v0.3.0

PhysioNet path: `https://physionet.org/files/mimic-iv-ecg/0.3.0/` ([physionet.org](https://physionet.org/content/mimic-iv-ecg/?utm_source=chatgpt.com))

```plaintext
mimic-iv-ecg/
└── 0.3.0/
    ├── 10045_0001/         # Record directory (subjectID_recordID)
    │   ├── 10045_0001.dat  # ECG data file (12-lead, 10s @500Hz)
    │   └── 10045_0001.hea  # ECG header file
    └── ...
```

Notes:

* Contains diagnostic 12‑lead ECG recordings matched to the MIMIC‑IV Clinical Database. ([physionet.org](https://physionet.org/content/mimic-iv-ecg/?utm_source=chatgpt.com))
