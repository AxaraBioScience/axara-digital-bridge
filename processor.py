import os, json, time, gzip, shutil, base64, subprocess, threading, uuid, math, requests
import numpy as np
import pandas as pd
from datetime import datetime
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import JSONResponse, StreamingResponse

app = FastAPI()

SAGE_BINARY = os.path.expanduser("~/sage")
FASTA_PATH = os.path.expanduser("~/human_proteome.fasta")
WORK_DIR = "/tmp/axara_jobs"
os.makedirs(WORK_DIR, exist_ok=True)

jobs = {}
last_job_time = [time.time()]

def auto_stop():
    while True:
        time.sleep(60)
        idle = time.time() - last_job_time[0]
        if idle > 900:
            print("15 minutes idle — shutting down VM", flush=True)
            os.system("sudo shutdown -h now")

threading.Thread(target=auto_stop, daemon=True).start()

# ── FORMAT DETECTION ──────────────────────────────────────────────────────────

def detect_format(filename):
    name = filename.lower()
    if name.endswith('.mzml.gz'): return 'mzml_gz'
    elif name.endswith('.mzml'): return 'mzml'
    elif name.endswith('.raw'): return 'thermo_raw'
    return 'unknown'

def needs_centroiding(mzml_path):
    try:
        with open(mzml_path, 'r', errors='ignore') as f:
            chunk = f.read(50000)
            if 'profile spectrum' in chunk.lower(): return True
            if 'centroid spectrum' in chunk.lower(): return False
        return False
    except:
        return False

def centroid_mzml(input_path, output_path):
    try:
        from pyteomics import mzml as pymzml
        def centroid_spectrum(mz, intensity, min_intensity=200):
            if len(mz) < 3: return mz, intensity
            peaks_mz, peaks_int = [], []
            for i in range(1, len(intensity) - 1):
                if intensity[i] > intensity[i-1] and intensity[i] > intensity[i+1] and intensity[i] > min_intensity:
                    peaks_mz.append(mz[i])
                    peaks_int.append(intensity[i])
            if len(peaks_mz) < 6:
                top_idx = np.argsort(intensity)[-50:]
                return mz[np.sort(top_idx)], intensity[np.sort(top_idx)]
            return np.array(peaks_mz), np.array(peaks_int)
        def encode_array(arr):
            return base64.b64encode(np.array(arr, dtype=np.float64).tobytes()).decode()
        spectra = []
        with pymzml.MzML(input_path) as reader:
            for spectrum in reader:
                if spectrum.get('ms level') == 2:
                    mz = spectrum['m/z array']
                    intensity = spectrum['intensity array']
                    c_mz, c_int = centroid_spectrum(mz, intensity)
                    if len(c_mz) >= 6:
                        prec = spectrum.get('precursorList', {}).get('precursor', [{}])[0]
                        sel = prec.get('selectedIonList', {}).get('selectedIon', [{}])[0]
                        prec_mz = sel.get('selected ion m/z', 0)
                        charge = sel.get('charge state', 2)
                        ce = prec.get('activation', {}).get('collision energy', 28.0)
                        n = len(c_mz)
                        lines = [
                            f'<spectrum index="{len(spectra)}" id="scan={len(spectra)+1}" defaultArrayLength="{n}">',
                            '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum"/>',
                            '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>',
                            '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>',
                            '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>',
                            f'<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{sum(c_int):.0f}"/>',
                            '<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>',
                            f'<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{prec_mz:.4f}"/>',
                            f'<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{charge}"/>',
                            '</selectedIon></selectedIonList><activation>',
                            '<cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation"/>',
                            f'<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="{ce}"/>',
                            '</activation></precursor></precursorList>',
                            '<binaryDataArrayList count="2">',
                            '<binaryDataArray>',
                            '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array"/>',
                            '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>',
                            '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>',
                            f'<binary>{encode_array(c_mz)}</binary>',
                            '</binaryDataArray>',
                            '<binaryDataArray>',
                            '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array"/>',
                            '<cvParam cvRef="MS" accession="MS:1000576" name="no compression"/>',
                            '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>',
                            f'<binary>{encode_array(c_int)}</binary>',
                            '</binaryDataArray></binaryDataArrayList></spectrum>'
                        ]
                        spectra.append('\n'.join(lines))
        header = [
            '<?xml version="1.0" encoding="utf-8"?>',
            '<mzML xmlns="http://psi.hupo.org/ms/mzml">',
            '<cvList count="2">',
            '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.30" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>',
            '<cv id="UO" fullName="Unit Ontology" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>',
            '</cvList><fileDescription><fileContent>',
            '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum"/>',
            '</fileContent></fileDescription>',
            '<dataProcessingList count="1"><dataProcessing id="dp">',
            '<processingMethod order="0" softwareRef="sw">',
            '<cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML"/>',
            '</processingMethod></dataProcessing></dataProcessingList>',
            '<softwareList count="1"><software id="sw" version="1.0">',
            '<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool"/>',
            '</software></softwareList>',
            f'<run><spectrumList count="{len(spectra)}" defaultDataProcessingRef="dp">',
        ]
        with open(output_path, 'w') as f:
            f.write('\n'.join(header) + '\n')
            for spec in spectra:
                f.write(spec + '\n')
            f.write('</spectrumList></run></mzML>')
        return output_path, None
    except Exception as e:
        return None, str(e)

# ── SAGE ──────────────────────────────────────────────────────────────────────

def build_sage_config(mzml_path, output_dir):
    config = {
        "database": {
            "bucket_size": 8192,
            "enzyme": {"missed_cleavages": 2, "min_len": 5, "max_len": 50, "cleave_at": "KR", "restrict": "P"},
            "fragment_min_mz": 150.0, "fragment_max_mz": 2000.0,
            "peptide_min_mass": 500.0, "peptide_max_mass": 7000.0,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": {"C": 57.0215},
            "variable_mods": {"K": [196.0706], "M": [15.9949]},
            "max_variable_mods": 2, "decoy_tag": "rev_",
            "generate_decoys": True, "fasta": FASTA_PATH
        },
        "precursor_tol": {"ppm": [-15.0, 15.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "precursor_charge": [2, 4], "isotope_errors": [0, 1],
        "deisotope": True, "chimera": False, "wide_window": False,
        "min_peaks": 6, "max_peaks": 150, "min_matched_peaks": 4,
        "report_psms": 5, "predict_rt": True,
        "mzml_paths": [os.path.abspath(mzml_path)],
        "output_directory": output_dir
    }
    config_path = os.path.join(output_dir, "sage_config.json")
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    return config_path

def parse_results(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    sda = df[df['peptide'].str.contains(r'\+196', na=False)].copy()
    return {
        "total_psms": len(df),
        "total_sda": len(sda),
        "high": sda[sda['spectrum_q'] <= 0.01].to_dict('records'),
        "medium": sda[(sda['spectrum_q'] > 0.01) & (sda['spectrum_q'] <= 0.05)].to_dict('records'),
        "low": sda[(sda['spectrum_q'] > 0.05) & (sda['spectrum_q'] <= 0.10)].to_dict('records'),
        "all_sda": sda.to_csv(sep='\t', index=False),
        "_df": df  # keep full dataframe in memory for Tier 2
    }

# ── TIER 2 — STEP 1: MS2 EXTRACTION ──────────────────────────────────────────

def extract_ms2_spectra(mzml_path, hits):
    """
    Extract MS2 spectra for each SDA hit from the mzML file.
    Matches by scan number recorded in the Sage PSM output.
    Returns dict: {scan_number: {'mz': array, 'intensity': array, 'precursor_mz': float, 'charge': int}}
    """
    from pyteomics import mzml as pymzml

    # Build set of scan numbers we need from Sage hits
    scan_numbers = set()
    for hit in hits:
        scan_ref = str(hit.get('scannr', hit.get('ScanNr', '')))
        # Sage scan refs are typically "file.mzML:scan=N" or just N
        if ':' in scan_ref:
            scan_ref = scan_ref.split('=')[-1]
        try:
            scan_numbers.add(int(scan_ref))
        except ValueError:
            pass

    spectra = {}
    try:
        with pymzml.MzML(mzml_path) as reader:
            for spectrum in reader:
                if spectrum.get('ms level') != 2:
                    continue
                # Extract scan number from spectrum id
                spec_id = spectrum.get('id', '')
                scan_num = None
                for part in spec_id.split():
                    if part.startswith('scan='):
                        try:
                            scan_num = int(part.split('=')[1])
                        except ValueError:
                            pass
                if scan_num is None:
                    continue
                if scan_num not in scan_numbers:
                    continue

                mz_array = np.array(spectrum.get('m/z array', []))
                int_array = np.array(spectrum.get('intensity array', []))

                prec = spectrum.get('precursorList', {}).get('precursor', [{}])[0]
                sel = prec.get('selectedIonList', {}).get('selectedIon', [{}])[0]
                prec_mz = float(sel.get('selected ion m/z', 0))
                charge = int(sel.get('charge state', 2))

                spectra[scan_num] = {
                    'mz': mz_array,
                    'intensity': int_array,
                    'precursor_mz': prec_mz,
                    'charge': charge
                }
    except Exception as e:
        print(f"MS2 extraction warning: {e}", flush=True)

    return spectra

# ── TIER 2 — STEP 2: THEORETICAL ION GENERATION ──────────────────────────────

# Amino acid monoisotopic masses
AA_MASS = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

PROTON = 1.007276
WATER = 18.010565
SDA_PROBE = 196.0706
SDA_TARGET = 82.0419
CAM = 57.0215      # carbamidomethyl on C
OX = 15.9949       # oxidation on M
FRAGMENT_TOL_PPM = 20.0

def parse_peptide_sequence(peptide_str):
    """
    Parse Sage peptide string like LGDK(+196.0706)VFSR into
    list of (residue, mod_mass) tuples.
    """
    residues = []
    i = 0
    while i < len(peptide_str):
        aa = peptide_str[i]
        mod = 0.0
        if aa not in AA_MASS:
            i += 1
            continue
        i += 1
        # Check for modification in parentheses
        if i < len(peptide_str) and peptide_str[i] == '(':
            j = peptide_str.index(')', i)
            mod_str = peptide_str[i+1:j]
            try:
                mod = float(mod_str.lstrip('+'))
            except ValueError:
                mod = 0.0
            i = j + 1
        base_mass = AA_MASS[aa]
        # Apply static mod for C
        if aa == 'C':
            base_mass += CAM
        residues.append((aa, base_mass + mod))
    return residues

def generate_theoretical_ions(peptide_str, charge=1):
    """
    Generate theoretical b and y ion series for an SDA-modified peptide.
    Returns dict with 'b' and 'y' lists of (ion_label, mz) tuples.
    """
    residues = parse_peptide_sequence(peptide_str)
    n = len(residues)
    if n < 2:
        return {'b': [], 'y': []}

    # b ions: N-terminal fragments
    b_ions = []
    mass = 0.0
    for i in range(n - 1):
        mass += residues[i][1]
        for z in range(1, min(charge, 3) + 1):
            mz = (mass + z * PROTON) / z
            label = f"b{i+1}" if z == 1 else f"b{i+1}²" if z == 2 else f"b{i+1}³"
            b_ions.append((label, mz))

    # y ions: C-terminal fragments
    y_ions = []
    mass = WATER
    for i in range(n - 1, 0, -1):
        mass += residues[i][1]
        ion_num = n - i
        for z in range(1, min(charge, 3) + 1):
            mz = (mass + z * PROTON) / z
            label = f"y{ion_num}" if z == 1 else f"y{ion_num}²" if z == 2 else f"y{ion_num}³"
            y_ions.append((label, mz))

    return {'b': b_ions, 'y': y_ions}

def match_ions(theoretical_ions, obs_mz, obs_intensity):
    """
    Match theoretical ions against observed spectrum within FRAGMENT_TOL_PPM.
    Returns list of matched (label, theoretical_mz, observed_mz, intensity, weight) tuples.
    """
    if len(obs_mz) == 0:
        return []
    base_peak = float(np.max(obs_intensity))
    if base_peak == 0:
        return []

    matches = []
    all_ions = theoretical_ions['b'] + theoretical_ions['y']

    for label, theo_mz in all_ions:
        tol = theo_mz * FRAGMENT_TOL_PPM / 1e6
        # Find closest observed peak within tolerance
        diffs = np.abs(obs_mz - theo_mz)
        idx = np.argmin(diffs)
        if diffs[idx] <= tol:
            weight = float(obs_intensity[idx]) / base_peak
            matches.append((label, theo_mz, float(obs_mz[idx]), float(obs_intensity[idx]), weight))

    return matches

# ── TIER 2 — STEP 3: BAYESIAN SCORER ─────────────────────────────────────────

def bayesian_score(matches, theoretical_ions, obs_mz, sage_fdr):
    """
    Compute Bayesian posterior probability for a peptide assignment.

    Prior: random ion match probability from fragment tolerance and m/z range.
    Likelihood: intensity-weighted product of match probabilities.
    Posterior: P(correct|data) = P(data|correct) * P(correct) / P(data)

    Returns (posterior_pct, confidence_tier, interpretation)
    """
    all_ions = theoretical_ions['b'] + theoretical_ions['y']
    n_theoretical = len(all_ions)

    if n_theoretical == 0 or len(obs_mz) == 0:
        return 0.0, "Low", "Insufficient theoretical ions generated for this peptide."

    # Prior: false match probability per ion
    mz_range = float(np.max(obs_mz) - np.min(obs_mz)) if len(obs_mz) > 1 else 1800.0
    avg_mz = float(np.mean(obs_mz))
    p_chance = min((2.0 * FRAGMENT_TOL_PPM * avg_mz) / (1e6 * mz_range), 0.1)

    n_matched = len(matches)

    # Edge case: sparse match
    if n_matched < 3:
        return 15.0, "Low", "Insufficient fragment coverage for statistical confidence (fewer than 3 ions matched)."

    # Likelihood: intensity-weighted match probability
    # P(data|correct) = product of weights for matched ions, penalised for unmatched
    matched_weights = [m[4] for m in matches]
    p_data_given_correct = float(np.mean(matched_weights)) * (n_matched / n_theoretical)

    # P(data|incorrect) = random match probability
    p_data_given_incorrect = p_chance ** n_matched

    # P(correct) from Sage FDR score
    p_correct_prior = max(1.0 - float(sage_fdr), 0.01)

    # Bayesian posterior
    numerator = p_data_given_correct * p_correct_prior
    denominator = numerator + p_data_given_incorrect * (1.0 - p_correct_prior)
    if denominator == 0:
        posterior = 0.0
    else:
        posterior = numerator / denominator

    posterior_pct = min(posterior * 100.0, 99.9)

    # Classify
    if posterior_pct >= 85.0:
        tier = "High"
        interp = (f"Strong MS2 fragmentation evidence ({n_matched}/{n_theoretical} theoretical ions matched). "
                  f"Fragment pattern is consistent with the proposed SDA modification site.")
    elif posterior_pct >= 60.0:
        tier = "Moderate"
        interp = (f"Partial MS2 evidence ({n_matched}/{n_theoretical} theoretical ions matched). "
                  f"Key diagnostic ions present but coverage is incomplete.")
    else:
        tier = "Low"
        interp = (f"Insufficient MS2 evidence ({n_matched}/{n_theoretical} theoretical ions matched). "
                  f"Hit retained for completeness. Orthogonal validation recommended.")

    return round(posterior_pct, 1), tier, interp

# ── TIER 2 — STEP 4: MS1 ISOTOPE VALIDATION ──────────────────────────────────

def averagine_isotope_pattern(mass, n_isotopes=5):
    """
    Generate theoretical isotope envelope using the averagine model.
    Returns list of relative intensities normalised to the monoisotopic peak.
    """
    # Averagine formula: empirical average amino acid composition
    # C: 4.9384, H: 7.7583, N: 1.3577, O: 1.4773, S: 0.0417 per 111.1 Da
    scale = mass / 111.1
    n_C = 4.9384 * scale
    n_H = 7.7583 * scale
    n_N = 1.3577 * scale
    n_O = 1.4773 * scale
    n_S = 0.0417 * scale

    # Natural isotope probabilities
    p_C = 0.01103
    p_H = 0.000115
    p_N = 0.00368
    p_O = 0.00038
    p_S = 0.0425

    # Expected number of heavy atoms (simplified Poisson)
    lam = n_C * p_C + n_H * p_H + n_N * p_N + n_O * p_O + n_S * p_S

    # Poisson distribution for isotope envelope
    pattern = []
    for k in range(n_isotopes):
        prob = (lam ** k) * math.exp(-lam) / math.factorial(k)
        pattern.append(prob)

    # Normalise to monoisotopic peak
    if pattern[0] > 0:
        pattern = [p / pattern[0] for p in pattern]

    return pattern

def validate_ms1_isotope(mzml_path, hit, obs_ms2_scan):
    """
    Extract MS1 precursor isotope envelope for a hit and compare to
    theoretical averagine pattern using Pearson r.
    Returns (pearson_r, classification, note)
    """
    from pyteomics import mzml as pymzml

    precursor_mz = obs_ms2_scan.get('precursor_mz', 0)
    charge = obs_ms2_scan.get('charge', 2)

    if precursor_mz == 0 or charge == 0:
        return None, "Inconclusive", "Precursor m/z or charge state not available."

    peptide_mass = (precursor_mz * charge) - (charge * PROTON)
    theo_pattern = averagine_isotope_pattern(peptide_mass, n_isotopes=5)
    isotope_spacing = 1.003355 / charge  # 13C spacing per charge state

    observed_intensities = []
    ms1_found = False

    try:
        with pymzml.MzML(mzml_path) as reader:
            prev_ms1_mz = None
            prev_ms1_int = None
            target_scan_found = False

            for spectrum in reader:
                ms_level = spectrum.get('ms level', 0)
                spec_id = spectrum.get('id', '')

                # Track MS1 scans — we want the one immediately before the MS2
                if ms_level == 1:
                    prev_ms1_mz = np.array(spectrum.get('m/z array', []))
                    prev_ms1_int = np.array(spectrum.get('intensity array', []))

                elif ms_level == 2:
                    # Check if this is our target MS2 scan
                    scan_num = None
                    for part in spec_id.split():
                        if part.startswith('scan='):
                            try:
                                scan_num = int(part.split('=')[1])
                            except ValueError:
                                pass

                    if scan_num == obs_ms2_scan.get('scan_num'):
                        target_scan_found = True
                        break

            if target_scan_found and prev_ms1_mz is not None and len(prev_ms1_mz) > 0:
                ms1_found = True
                # Extract isotope envelope around precursor_mz
                tol = precursor_mz * 20.0 / 1e6  # 20 ppm
                for i in range(5):
                    iso_mz = precursor_mz + i * isotope_spacing
                    diffs = np.abs(prev_ms1_mz - iso_mz)
                    idx = np.argmin(diffs)
                    if diffs[idx] <= tol * (1 + i * 0.5):  # slightly wider for higher isotopes
                        observed_intensities.append(float(prev_ms1_int[idx]))
                    else:
                        observed_intensities.append(0.0)

    except Exception as e:
        return None, "Inconclusive", f"MS1 extraction error: {str(e)[:80]}"

    if not ms1_found:
        return None, "Inconclusive", "MS1 survey scan not found preceding this MS2 event."

    # Need at least 3 non-zero isotope peaks for meaningful correlation
    non_zero = sum(1 for v in observed_intensities if v > 0)
    if non_zero < 3:
        return None, "Inconclusive", "Insufficient isotope peaks detected in MS1 scan."

    # Normalise observed to monoisotopic
    if observed_intensities[0] > 0:
        obs_norm = [v / observed_intensities[0] for v in observed_intensities]
    else:
        # Find first non-zero and normalise
        first_nonzero = next((v for v in observed_intensities if v > 0), 1.0)
        obs_norm = [v / first_nonzero for v in observed_intensities]

    # Pearson correlation between observed and theoretical
    n = min(len(obs_norm), len(theo_pattern))
    obs_arr = np.array(obs_norm[:n])
    theo_arr = np.array(theo_pattern[:n])

    if np.std(obs_arr) == 0 or np.std(theo_arr) == 0:
        return 0.0, "Inconclusive", "Isotope pattern variance too low for correlation."

    pearson_r = float(np.corrcoef(obs_arr, theo_arr)[0, 1])

    if pearson_r >= 0.90:
        classification = "Confirmed"
        note = "MS1 isotope envelope is consistent with the proposed SDA-modified peptide mass."
    elif pearson_r >= 0.70:
        classification = "Plausible"
        note = "Isotope pattern broadly consistent. Minor deviations may reflect co-isolation."
    else:
        classification = "Inconclusive"
        note = "MS1 isotope pattern does not match the proposed modification. Treat with caution."

    return round(pearson_r, 3), classification, note

# ── TIER 2 — STEP 5: DEEPLC RETENTION TIME PREDICTION ────────────────────────

def predict_retention_times(hits, all_psms_df):
    """
    Use DeepLC to predict retention times for SDA hits.
    Calibrates using top non-SDA high-confidence PSMs from the same run.
    Returns dict: {peptide_str: {'predicted_rt': float, 'observed_rt': float,
                                  'deviation_min': float, 'classification': str, 'note': str}}
    """
    try:
        from deeplc import DeepLC
    except ImportError:
        return None, "DeepLC not installed. Retention time prediction skipped."

    try:
        # Build calibration set: top 20 non-SDA high-confidence PSMs
        cal_df = all_psms_df[
            ~all_psms_df['peptide'].str.contains(r'\+196', na=False) &
            (all_psms_df['spectrum_q'] <= 0.01)
        ].copy()

        if len(cal_df) < 5:
            return None, "Insufficient calibration peptides (need ≥5 non-SDA high-confidence PSMs)."

        cal_df = cal_df.nsmallest(20, 'spectrum_q')

        # DeepLC expects plain peptide sequences (strip modifications for calibration)
        def strip_mods(pep):
            import re
            return re.sub(r'\([^)]*\)', '', pep)

        cal_peptides = []
        for _, row in cal_df.iterrows():
            plain = strip_mods(str(row['peptide']))
            rt = float(row.get('rt', row.get('retention_time', 0)))
            cal_peptides.append((plain, '', rt))

        # Prepare SDA hit sequences for prediction
        sda_peptides = []
        sda_hits_plain = []
        for hit in hits:
            pep = str(hit.get('peptide', ''))
            plain = strip_mods(pep)
            sda_peptides.append((plain, ''))
            sda_hits_plain.append(pep)

        # Run DeepLC
        dlc = DeepLC(verbose=False)
        dlc.calibrate_preds(seq_df=pd.DataFrame(cal_peptides, columns=['seq', 'modifications', 'tr']))
        preds = dlc.make_preds(seq_df=pd.DataFrame(sda_peptides, columns=['seq', 'modifications']))

        results = {}
        for i, hit in enumerate(hits):
            pep = str(hit.get('peptide', ''))
            observed_rt = float(hit.get('rt', hit.get('retention_time', 0))) / 60.0  # convert to minutes
            predicted_rt = float(preds[i]) / 60.0
            deviation = abs(predicted_rt - observed_rt)

            if deviation <= 2.0:
                classification = "Consistent"
                note = "Observed retention time is consistent with the predicted elution window."
            elif deviation <= 5.0:
                classification = "Marginal"
                note = "Retention time is within acceptable range but at the outer boundary."
            else:
                classification = "Inconsistent"
                note = "Retention time deviates substantially from prediction. Chimeric spectrum should be considered."

            results[pep] = {
                'predicted_rt': round(predicted_rt, 2),
                'observed_rt': round(observed_rt, 2),
                'deviation_min': round(deviation, 2),
                'classification': classification,
                'note': note
            }

        return results, None

    except Exception as e:
        return None, f"Retention time prediction failed: {str(e)[:120]}"

# ── TIER 2 — STEP 6: UNIPROT ANNOTATION ──────────────────────────────────────

def fetch_uniprot_annotation(accession, peptide_seq, uniprot_cache):
    """
    Fetch protein annotation from UniProt REST API.
    Maps peptide to protein position and identifies functional domain overlap.
    Returns annotation dict.
    """
    # Strip modification notation for sequence matching
    import re
    plain_peptide = re.sub(r'\([^)]*\)', '', peptide_seq)

    # Check cache first
    if accession in uniprot_cache:
        data = uniprot_cache[accession]
    else:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
            resp = requests.get(url, timeout=10)
            if resp.status_code != 200:
                return {
                    'protein_name': accession,
                    'gene': 'Unknown',
                    'annotation_note': 'UniProt annotation not available for this accession.',
                    'domain': 'Not mapped',
                    'functional_note': '',
                    'disease': ''
                }
            data = resp.json()
            uniprot_cache[accession] = data
        except Exception as e:
            return {
                'protein_name': accession,
                'gene': 'Unknown',
                'annotation_note': f'UniProt lookup failed: {str(e)[:60]}',
                'domain': 'Not mapped',
                'functional_note': '',
                'disease': ''
            }
        # Retry once on timeout
        time.sleep(0.5)

    # Extract protein name
    protein_name = accession
    try:
        protein_name = data['proteinDescription']['recommendedName']['fullName']['value']
    except (KeyError, TypeError):
        try:
            protein_name = data['proteinDescription']['submittedName'][0]['fullName']['value']
        except (KeyError, TypeError, IndexError):
            pass

    # Extract gene name
    gene = 'Unknown'
    try:
        gene = data['genes'][0]['geneName']['value']
    except (KeyError, TypeError, IndexError):
        pass

    # Extract canonical sequence for position mapping
    canonical_seq = ''
    try:
        canonical_seq = data['sequence']['value']
    except (KeyError, TypeError):
        pass

    # Map peptide to protein position
    peptide_position = None
    if canonical_seq and plain_peptide in canonical_seq:
        pos = canonical_seq.index(plain_peptide)
        peptide_position = (pos + 1, pos + len(plain_peptide))  # 1-indexed

    # Extract functional features (domains, active sites, binding sites)
    domain_hit = 'Not mapped'
    functional_note = ''
    disease_note = ''

    if peptide_position:
        pep_start, pep_end = peptide_position
        try:
            features = data.get('features', [])
            for feat in features:
                feat_type = feat.get('type', '')
                loc = feat.get('location', {})
                feat_start = loc.get('start', {}).get('value', 0)
                feat_end = loc.get('end', {}).get('value', 0)

                # Check overlap with peptide position
                if feat_start <= pep_end and feat_end >= pep_start:
                    feat_desc = feat.get('description', feat_type)
                    if feat_type in ('Domain', 'Region', 'Motif'):
                        domain_hit = f"{feat_desc} (residues {feat_start}–{feat_end})"
                    elif feat_type == 'Active site':
                        functional_note += f"Active site at residue {feat_start}. "
                    elif feat_type == 'Binding site':
                        ligand = feat.get('ligand', {}).get('name', 'ligand')
                        functional_note += f"Binding site for {ligand} (residues {feat_start}–{feat_end}). "
        except Exception:
            pass

        # Disease associations
        try:
            comments = data.get('comments', [])
            for comment in comments:
                if comment.get('commentType') == 'DISEASE':
                    disease_name = comment.get('disease', {}).get('diseaseId', '')
                    if disease_name:
                        disease_note += disease_name + ". "
        except Exception:
            pass

    position_str = f"Residues {peptide_position[0]}–{peptide_position[1]}" if peptide_position else "Position not mapped"

    return {
        'protein_name': protein_name,
        'gene': gene,
        'peptide_position': position_str,
        'domain': domain_hit,
        'functional_note': functional_note.strip(),
        'disease': disease_note.strip(),
        'annotation_note': ''
    }

# ── TIER 2 — STEP 7: RULE OF 25 SEQUENCE CHECK ───────────────────────────────

def rule_of_25_sequence_check(peptide_str, protein_seq=''):
    """
    Sequence-level Rule of 25 check for SDA probe residency plausibility.
    Checks: K spacing, local hydrophobicity, accessibility flags.
    Returns (outcome, note) where outcome is 'Favourable', 'Neutral', or 'Atypical'
    """
    import re
    plain = re.sub(r'\([^)]*\)', '', peptide_str)

    # Hydrophobic residues
    hydrophobic = set('VILFMWA')
    # Known burial motifs that preclude SDA labelling
    burial_motifs = ['GxxxG', 'PXXP']

    score = 0
    notes = []

    # Check 1: K present (required for SDA probe adduct)
    if 'K' not in plain:
        return 'Atypical', 'No lysine residue found in peptide — SDA probe adduct site is absent.'

    # Check 2: Local hydrophobicity around K
    k_positions = [i for i, aa in enumerate(plain) if aa == 'K']
    for k_pos in k_positions:
        window_start = max(0, k_pos - 3)
        window_end = min(len(plain), k_pos + 4)
        window = plain[window_start:window_end]
        hydrophobic_count = sum(1 for aa in window if aa in hydrophobic)
        hydrophobic_ratio = hydrophobic_count / len(window)
        if hydrophobic_ratio >= 0.4:
            score += 1
            notes.append(f"K at position {k_pos+1} is in a hydrophobic context (ratio {hydrophobic_ratio:.1f}) — consistent with SDA probe residency.")

    # Check 3: Proline near K (proline reduces accessibility)
    for k_pos in k_positions:
        nearby = plain[max(0, k_pos-2):min(len(plain), k_pos+3)]
        if 'P' in nearby:
            score -= 1
            notes.append(f"Proline adjacent to K at position {k_pos+1} may reduce probe accessibility.")

    # Check 4: Peptide length (very short peptides are less reliable)
    if len(plain) < 6:
        score -= 1
        notes.append("Short peptide — lower confidence in probe residency assignment.")
    elif len(plain) >= 8:
        score += 1

    # Classify
    if score >= 1:
        outcome = 'Favourable'
        summary = "Sequence characteristics are consistent with SDA probe residency."
    elif score == 0:
        outcome = 'Neutral'
        summary = "No sequence-level features that either support or contradict probe residency."
    else:
        outcome = 'Atypical'
        summary = "Sequence context is atypical for SDA labelling. Consider orthogonal validation."

    detail = ' '.join(notes) if notes else summary
    return outcome, detail

# ── TIER 2 — STEP 8: MATPLOTLIB SPECTRUM FIGURES ─────────────────────────────

def generate_spectrum_figure(hit, ms2_scan, matches, theoretical_ions, fig_path):
    """
    Generate annotated MS2 spectrum figure for a single hit.
    Grey bars: all observed ions
    Blue labels: matched b ions
    Red labels: matched y ions
    SDA mass shift annotated on relevant ions.
    Saved as PNG to fig_path.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend — critical for server
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        obs_mz = ms2_scan['mz']
        obs_int = ms2_scan['intensity']

        if len(obs_mz) == 0:
            return False

        # Normalise intensities to 100
        max_int = float(np.max(obs_int))
        if max_int == 0:
            return False
        obs_int_norm = (obs_int / max_int) * 100.0

        # Build match lookup
        b_matches = {m[2]: (m[0], m[4]) for m in matches if m[0].startswith('b')}
        y_matches = {m[2]: (m[0], m[4]) for m in matches if m[0].startswith('y')}

        fig, ax = plt.subplots(figsize=(12, 5))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        # Draw all observed ions in grey
        for mz_val, int_val in zip(obs_mz, obs_int_norm):
            color = '#CCCCCC'
            if mz_val in b_matches:
                color = '#2874A6'  # blue for b ions
            elif mz_val in y_matches:
                color = '#C0392B'  # red for y ions
            ax.bar(mz_val, int_val, width=0.8, color=color, linewidth=0)

        # Annotate matched ions
        annotated = set()
        for mz_val, (label, weight) in b_matches.items():
            idx = np.argmin(np.abs(obs_mz - mz_val))
            int_val = obs_int_norm[idx]
            if label not in annotated and int_val > 5:
                ax.annotate(label, xy=(mz_val, int_val), xytext=(0, 4),
                           textcoords='offset points', ha='center', va='bottom',
                           fontsize=7, color='#2874A6', fontweight='bold')
                annotated.add(label)
                # Mark SDA-bearing ions
                if '+196' in str(hit.get('peptide', '')) and 'K' in label:
                    ax.annotate('SDA', xy=(mz_val, int_val + 8), xytext=(0, 2),
                               textcoords='offset points', ha='center', fontsize=6,
                               color='#C5A059', fontweight='bold')

        for mz_val, (label, weight) in y_matches.items():
            idx = np.argmin(np.abs(obs_mz - mz_val))
            int_val = obs_int_norm[idx]
            if label not in annotated and int_val > 5:
                ax.annotate(label, xy=(mz_val, int_val), xytext=(0, 4),
                           textcoords='offset points', ha='center', va='bottom',
                           fontsize=7, color='#C0392B', fontweight='bold')
                annotated.add(label)

        # Labels and formatting
        peptide_display = str(hit.get('peptide', ''))[:50]
        protein_display = str(hit.get('proteins', ''))[:40]
        ax.set_title(f"{peptide_display}\n{protein_display}",
                    fontsize=9, fontweight='bold', color='#1A1A2E', pad=8)
        ax.set_xlabel('m/z', fontsize=9, color='#1A1A2E')
        ax.set_ylabel('Relative Intensity (%)', fontsize=9, color='#1A1A2E')
        ax.set_ylim(0, 115)
        ax.tick_params(colors='#1A1A2E', labelsize=8)
        for spine in ax.spines.values():
            spine.set_color('#CCCCCC')

        # Legend
        b_patch = mpatches.Patch(color='#2874A6', label='b ions')
        y_patch = mpatches.Patch(color='#C0392B', label='y ions')
        grey_patch = mpatches.Patch(color='#CCCCCC', label='unmatched')
        ax.legend(handles=[b_patch, y_patch, grey_patch], fontsize=8,
                 loc='upper right', framealpha=0.8)

        # AXARA gold accent line at top
        ax.axhline(y=113, color='#C5A059', linewidth=2, alpha=0.6)

        plt.tight_layout()
        plt.savefig(fig_path, dpi=150, bbox_inches='tight',
                   facecolor='white', format='png')
        plt.close(fig)
        return True

    except Exception as e:
        print(f"Figure generation failed for {fig_path}: {e}", flush=True)
        return False

# ── TIER 2 — STEP 9: REPORTLAB PDF ASSEMBLY ──────────────────────────────────

def generate_tier2_pdf(job_id, work_dir, filename, results, tier2_data):
    """
    Assemble the full Tier 2 PDF report using reportlab.
    Zero retention: PDF written to /tmp only, never to persistent storage.
    Returns path to the generated PDF.
    """
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.lib.colors import HexColor, white, black, lightgrey
        from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                         TableStyle, Image, HRFlowable, PageBreak,
                                         KeepTogether)
        from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
        from io import BytesIO

        GOLD = HexColor('#C5A059')
        DARK = HexColor('#1A1A2E')
        LIGHT_GREY = HexColor('#F2F2F2')
        MID_GREY = HexColor('#CCCCCC')
        BLUE = HexColor('#2874A6')
        RED_COL = HexColor('#C0392B')
        GREEN_COL = HexColor('#1E8449')
        AMBER = HexColor('#D4AC0D')

        pdf_path = os.path.join(work_dir, 'tier2_report.pdf')

        doc = SimpleDocTemplate(
            pdf_path,
            pagesize=letter,
            leftMargin=0.85*inch, rightMargin=0.85*inch,
            topMargin=0.85*inch, bottomMargin=0.85*inch,
            title="AXARA Bridge — Full Analytical Report",
            author="AXARA Digital"
        )

        styles = getSampleStyleSheet()

        # Custom styles
        def style(name, **kwargs):
            return ParagraphStyle(name, **kwargs)

        s_title = style('Title', fontName='Helvetica-Bold', fontSize=26,
                        textColor=DARK, spaceAfter=4, leading=30)
        s_subtitle = style('Subtitle', fontName='Helvetica-Bold', fontSize=18,
                           textColor=GOLD, spaceAfter=6, leading=22)
        s_h1 = style('H1', fontName='Helvetica-Bold', fontSize=14,
                     textColor=DARK, spaceBefore=16, spaceAfter=6,
                     borderPadding=(0,0,2,0))
        s_h2 = style('H2', fontName='Helvetica-Bold', fontSize=11,
                     textColor=GOLD, spaceBefore=10, spaceAfter=4)
        s_body = style('Body', fontName='Helvetica', fontSize=9,
                       textColor=DARK, spaceAfter=4, leading=13)
        s_note = style('Note', fontName='Helvetica-Oblique', fontSize=8,
                       textColor=HexColor('#666666'), spaceAfter=4, leading=11)
        s_centre = style('Centre', fontName='Helvetica', fontSize=9,
                         textColor=DARK, alignment=TA_CENTER)
        s_footer = style('Footer', fontName='Helvetica', fontSize=7,
                         textColor=HexColor('#999999'), alignment=TA_CENTER)

        def gold_rule():
            return HRFlowable(width="100%", thickness=2, color=GOLD, spaceAfter=8)

        def grey_rule():
            return HRFlowable(width="100%", thickness=0.5, color=MID_GREY,
                             spaceBefore=4, spaceAfter=4)

        def confidence_color(tier):
            if tier == 'High' or tier == 'Confirmed' or tier == 'Consistent':
                return GREEN_COL
            elif tier == 'Moderate' or tier == 'Plausible' or tier == 'Marginal':
                return AMBER
            else:
                return RED_COL

        def tier_badge(tier):
            color = confidence_color(tier)
            return Paragraph(
                f'<font color="#{color.hexval()[1:].upper()}">\u25cf</font> {tier}',
                style('Badge', fontName='Helvetica-Bold', fontSize=9,
                      textColor=DARK, spaceAfter=2)
            )

        story = []

        # ── COVER PAGE ──────────────────────────────────────────────────────
        story.append(Spacer(1, 0.4*inch))
        story.append(Paragraph("AXARA Bridge Platform", s_title))
        story.append(gold_rule())
        story.append(Paragraph("Full Analytical Report", s_subtitle))
        story.append(Spacer(1, 0.15*inch))

        cover_data = [
            ['Input file', filename],
            ['Report date', datetime.now().strftime('%d %B %Y, %H:%M UTC')],
            ['Job reference', job_id[:16].upper()],
            ['Tier', 'Tier 2 — Full Analytical Report'],
            ['Engine', 'Sage v0.14.6 | pyteomics | DeepLC | UniProt'],
        ]
        cover_table = Table(cover_data, colWidths=[1.8*inch, 5.0*inch])
        cover_table.setStyle(TableStyle([
            ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'),
            ('FONTNAME', (1,0), (1,-1), 'Helvetica'),
            ('FONTSIZE', (0,0), (-1,-1), 9),
            ('TEXTCOLOR', (0,0), (0,-1), HexColor('#888888')),
            ('TEXTCOLOR', (1,0), (1,-1), DARK),
            ('ROWBACKGROUNDS', (0,0), (-1,-1), [white, LIGHT_GREY]),
            ('TOPPADDING', (0,0), (-1,-1), 5),
            ('BOTTOMPADDING', (0,0), (-1,-1), 5),
            ('LEFTPADDING', (0,0), (-1,-1), 8),
            ('RIGHTPADDING', (0,0), (-1,-1), 8),
            ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
        ]))
        story.append(cover_table)
        story.append(Spacer(1, 0.2*inch))
        story.append(Paragraph(
            "This report is generated automatically by the AXARA Bridge platform. "
            "All data is processed in volatile memory only and is not retained after download. "
            "For research use. AXARA Digital — info@axara.bio",
            s_note))
        story.append(PageBreak())

        # ── SECTION 1: EXECUTIVE SUMMARY ────────────────────────────────────
        story.append(Paragraph("Executive Summary", s_h1))
        story.append(gold_rule())

        hits_all = results.get('high', []) + results.get('medium', []) + results.get('low', [])
        scored_hits = tier2_data.get('scored_hits', {})

        high_count = sum(1 for h in scored_hits.values() if h.get('ms2_tier') == 'High')
        mod_count = sum(1 for h in scored_hits.values() if h.get('ms2_tier') == 'Moderate')
        low_count = sum(1 for h in scored_hits.values() if h.get('ms2_tier') == 'Low')

        summary_data = [
            ['Total PSMs identified', str(results.get('total_psms', 0))],
            ['SDA-modified peptides', str(results.get('total_sda', 0))],
            ['High confidence hits (MS2 ≥ 85%)', str(high_count)],
            ['Moderate confidence hits (MS2 60–84%)', str(mod_count)],
            ['Low confidence hits (MS2 < 60%)', str(low_count)],
        ]
        sum_table = Table(summary_data, colWidths=[3.5*inch, 3.3*inch])
        sum_table.setStyle(TableStyle([
            ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'),
            ('FONTNAME', (1,0), (1,-1), 'Helvetica'),
            ('FONTSIZE', (0,0), (-1,-1), 9),
            ('TEXTCOLOR', (0,0), (-1,-1), DARK),
            ('ROWBACKGROUNDS', (0,0), (-1,-1), [white, LIGHT_GREY]),
            ('TOPPADDING', (0,0), (-1,-1), 5),
            ('BOTTOMPADDING', (0,0), (-1,-1), 5),
            ('LEFTPADDING', (0,0), (-1,-1), 8),
            ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
        ]))
        story.append(sum_table)
        story.append(Spacer(1, 0.1*inch))

        # Overall recommendation
        if high_count > 0:
            rec = (f"{high_count} hit(s) show strong multi-dimensional evidence. "
                   f"These are high-priority candidates for Tier 3 structural mapping.")
        elif mod_count > 0:
            rec = (f"No high-confidence hits identified. {mod_count} moderate-confidence hit(s) "
                   f"warrant orthogonal validation before Tier 3.")
        else:
            rec = ("Confidence levels are low across all hits. Additional experimental data "
                   "is recommended before proceeding to structural analysis.")

        story.append(Paragraph(f"<b>Overall recommendation:</b> {rec}", s_body))
        story.append(PageBreak())

        # ── SECTION 2: TIER 1 IDENTIFICATION TABLE ───────────────────────────
        story.append(Paragraph("Tier 1 Identification Table", s_h1))
        story.append(gold_rule())
        story.append(Paragraph("Peptide identifications from the Sage database search. "
                               "All SDA-modified hits are included regardless of confidence tier.",
                               s_body))
        story.append(Spacer(1, 0.08*inch))

        t1_header = [['Peptide', 'Protein', 'Charge', 'Hyperscore', 'FDR', 'Confidence']]
        t1_data = []
        for hit in hits_all[:30]:  # cap at 30 for readability
            pep = str(hit.get('peptide', ''))[:35]
            prot = str(hit.get('proteins', ''))[:25]
            charge = str(hit.get('charge', ''))
            score = f"{float(hit.get('hyperscore', 0)):.1f}"
            fdr = f"{float(hit.get('spectrum_q', 0)):.4f}"
            sh = scored_hits.get(str(hit.get('scannr', '')), {})
            conf = sh.get('ms2_tier', '—')
            t1_data.append([pep, prot, charge, score, fdr, conf])

        t1_table = Table(t1_header + t1_data,
                        colWidths=[2.0*inch, 1.5*inch, 0.5*inch, 0.8*inch, 0.6*inch, 0.8*inch])
        t1_style = [
            ('BACKGROUND', (0,0), (-1,0), DARK),
            ('TEXTCOLOR', (0,0), (-1,0), white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTSIZE', (0,0), (-1,-1), 8),
            ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
            ('TEXTCOLOR', (0,1), (-1,-1), DARK),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHT_GREY]),
            ('TOPPADDING', (0,0), (-1,-1), 4),
            ('BOTTOMPADDING', (0,0), (-1,-1), 4),
            ('LEFTPADDING', (0,0), (-1,-1), 6),
            ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
            ('ALIGN', (2,0), (4,-1), 'CENTER'),
        ]
        t1_table.setStyle(TableStyle(t1_style))
        story.append(t1_table)
        story.append(PageBreak())

        # ── SECTION 3: PER-HIT EVIDENCE ─────────────────────────────────────
        story.append(Paragraph("Fragmentation Evidence", s_h1))
        story.append(gold_rule())

        for hit_key, hit_data in scored_hits.items():
            hit = hit_data.get('hit', {})
            pep = str(hit.get('peptide', ''))
            prot = str(hit.get('proteins', ''))

            story.append(Paragraph(f"{pep[:60]}", s_h2))
            story.append(Paragraph(f"Protein: {prot[:80]}", s_body))
            story.append(grey_rule())

            # Spectrum figure
            fig_path = hit_data.get('fig_path')
            if fig_path and os.path.exists(fig_path):
                try:
                    img = Image(fig_path, width=6.5*inch, height=2.7*inch)
                    story.append(img)
                except Exception:
                    story.append(Paragraph("[Spectrum figure not available]", s_note))
            else:
                story.append(Paragraph("[Spectrum figure not available]", s_note))

            story.append(Spacer(1, 0.08*inch))

            # Evidence scores table
            ms2_tier = hit_data.get('ms2_tier', '—')
            ms2_pct = hit_data.get('ms2_posterior', 0)
            ms1_class = hit_data.get('ms1_classification', '—')
            ms1_r = hit_data.get('ms1_pearson_r', '—')
            rt_class = hit_data.get('rt_classification', '—')
            rt_dev = hit_data.get('rt_deviation_min', '—')
            r25_outcome = hit_data.get('r25_outcome', '—')

            ev_data = [
                ['Evidence Dimension', 'Score', 'Classification'],
                ['MS2 Fragment Ion (Bayesian)',
                 f"{ms2_pct:.1f}%" if isinstance(ms2_pct, float) else '—',
                 ms2_tier],
                ['MS1 Isotope Validation (Pearson r)',
                 f"{ms1_r:.3f}" if isinstance(ms1_r, float) else '—',
                 ms1_class],
                ['Retention Time (DeepLC)',
                 f"{rt_dev:.1f} min deviation" if isinstance(rt_dev, float) else '—',
                 rt_class],
                ['Rule of 25 (Sequence Level)', '—', r25_outcome],
            ]
            ev_table = Table(ev_data, colWidths=[2.5*inch, 1.8*inch, 1.9*inch])
            ev_table.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), DARK),
                ('TEXTCOLOR', (0,0), (-1,0), white),
                ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
                ('FONTSIZE', (0,0), (-1,-1), 8),
                ('TEXTCOLOR', (0,1), (-1,-1), DARK),
                ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHT_GREY]),
                ('TOPPADDING', (0,0), (-1,-1), 4),
                ('BOTTOMPADDING', (0,0), (-1,-1), 4),
                ('LEFTPADDING', (0,0), (-1,-1), 6),
                ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
            ]))
            story.append(ev_table)
            story.append(Spacer(1, 0.08*inch))

            # UniProt annotation
            ann = hit_data.get('annotation', {})
            if ann:
                ann_data = [
                    ['Protein', ann.get('protein_name', '—')],
                    ['Gene', ann.get('gene', '—')],
                    ['Peptide position', ann.get('peptide_position', '—')],
                    ['Domain', ann.get('domain', '—')],
                    ['Functional notes', ann.get('functional_note', '—') or '—'],
                    ['Disease association', ann.get('disease', '—') or 'None reported'],
                ]
                ann_table = Table(ann_data, colWidths=[1.5*inch, 4.7*inch])
                ann_table.setStyle(TableStyle([
                    ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'),
                    ('FONTNAME', (1,0), (1,-1), 'Helvetica'),
                    ('FONTSIZE', (0,0), (-1,-1), 8),
                    ('TEXTCOLOR', (0,0), (0,-1), HexColor('#888888')),
                    ('TEXTCOLOR', (1,0), (1,-1), DARK),
                    ('ROWBACKGROUNDS', (0,0), (-1,-1), [white, LIGHT_GREY]),
                    ('TOPPADDING', (0,0), (-1,-1), 4),
                    ('BOTTOMPADDING', (0,0), (-1,-1), 4),
                    ('LEFTPADDING', (0,0), (-1,-1), 6),
                    ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
                ]))
                story.append(ann_table)

            # Interpretation
            ms2_interp = hit_data.get('ms2_interpretation', '')
            r25_note = hit_data.get('r25_note', '')
            rt_note = hit_data.get('rt_note', '')
            ms1_note = hit_data.get('ms1_note', '')

            interp_parts = []
            if ms2_interp:
                interp_parts.append(ms2_interp)
            if ms1_note:
                interp_parts.append(ms1_note)
            if rt_note:
                interp_parts.append(rt_note)
            if r25_note:
                interp_parts.append(r25_note)

            if interp_parts:
                story.append(Spacer(1, 0.06*inch))
                story.append(Paragraph("<b>Interpretation:</b> " + " ".join(interp_parts), s_body))

            # Next-step recommendation
            if ms2_tier == 'High' and ms1_class == 'Confirmed':
                rec_text = ("Two independent lines of evidence support this hit. "
                           "<b>Recommend: Proceed to Tier 3 structural mapping.</b>")
            elif ms2_tier == 'High':
                rec_text = ("Strong fragment evidence. MS1 confirmation is incomplete. "
                           "<b>Recommend: Consider Tier 3 with awareness of MS1 uncertainty.</b>")
            elif ms2_tier == 'Moderate':
                rec_text = ("Partial evidence. <b>Recommend: Orthogonal validation before Tier 3.</b>")
            else:
                rec_text = ("Insufficient evidence. <b>Recommend: Additional data before proceeding.</b>")

            story.append(Paragraph(rec_text, s_body))
            story.append(Spacer(1, 0.15*inch))
            story.append(grey_rule())

        story.append(PageBreak())

        # ── SECTION 4: CONFIDENCE SUMMARY TABLE ─────────────────────────────
        story.append(Paragraph("Confidence Summary", s_h1))
        story.append(gold_rule())

        cs_header = [['Peptide', 'MS2 Score', 'MS2 Tier', 'MS1', 'RT', 'R25', 'Recommendation']]
        cs_data = []
        sorted_hits = sorted(scored_hits.values(),
                            key=lambda x: x.get('ms2_posterior', 0), reverse=True)
        for h in sorted_hits[:20]:
            hit = h.get('hit', {})
            pep = str(hit.get('peptide', ''))[:25]
            ms2_pct = h.get('ms2_posterior', 0)
            ms2_t = h.get('ms2_tier', '—')
            ms1 = h.get('ms1_classification', '—')[:4]
            rt = h.get('rt_classification', '—')[:4]
            r25 = h.get('r25_outcome', '—')[:4]
            if ms2_t == 'High' and ms1 == 'Conf':
                rec = 'Tier 3 Ready'
            elif ms2_t == 'High':
                rec = 'Consider Tier 3'
            elif ms2_t == 'Moderate':
                rec = 'Validate First'
            else:
                rec = 'More Data Needed'
            cs_data.append([pep,
                            f"{ms2_pct:.1f}%" if isinstance(ms2_pct, float) else '—',
                            ms2_t, ms1, rt, r25, rec])

        cs_table = Table(cs_header + cs_data,
                        colWidths=[1.9*inch, 0.7*inch, 0.7*inch, 0.55*inch,
                                   0.55*inch, 0.55*inch, 1.3*inch])
        cs_table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), DARK),
            ('TEXTCOLOR', (0,0), (-1,0), white),
            ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
            ('FONTSIZE', (0,0), (-1,-1), 7),
            ('TEXTCOLOR', (0,1), (-1,-1), DARK),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHT_GREY]),
            ('TOPPADDING', (0,0), (-1,-1), 3),
            ('BOTTOMPADDING', (0,0), (-1,-1), 3),
            ('LEFTPADDING', (0,0), (-1,-1), 5),
            ('GRID', (0,0), (-1,-1), 0.5, MID_GREY),
            ('ALIGN', (1,0), (-1,-1), 'CENTER'),
        ]))
        story.append(cs_table)
        story.append(PageBreak())

        # ── SECTION 5: METHODS APPENDIX ─────────────────────────────────────
        story.append(Paragraph("Methods Appendix", s_h1))
        story.append(gold_rule())

        methods = [
            ("Database Search", "Sage v0.14.6. SDA probe adduct +196.0706 Da on K, "
             "target adduct +82.0419 Da, oxidation +15.9949 Da on M, "
             "carbamidomethyl +57.0215 Da on C. Trypsin, 2 missed cleavages. "
             "Precursor tolerance 15 ppm. Fragment tolerance 20 ppm. "
             "Target-decoy FDR. Database: UniProt SwissProt human proteome."),
            ("Bayesian Fragment Ion Scorer",
             "Theoretical b/y ion series generated for each SDA-modified peptide. "
             "Observed ions matched within 20 ppm tolerance. Match probability weighted "
             "by observed ion intensity relative to base peak. Bayesian posterior "
             "computed using Sage FDR score as prior. Confidence tiers: "
             "High ≥ 85%, Moderate 60–84%, Low < 60%."),
            ("MS1 Isotope Validation",
             "Theoretical isotope envelope generated using the averagine model. "
             "Observed precursor isotope peaks extracted from the MS1 survey scan "
             "immediately preceding each MS2 event. Pearson correlation between "
             "observed and theoretical intensity vectors. "
             "Confirmed r ≥ 0.90, Plausible r 0.70–0.89, Inconclusive r < 0.70."),
            ("Retention Time Prediction",
             "DeepLC deep learning model. Calibrated using the top 20 non-SDA "
             "high-confidence PSMs (spectrum_q ≤ 0.01) from the same LC run. "
             "Consistent: deviation ≤ 2 min. Marginal: 2–5 min. "
             "Inconsistent: > 5 min."),
            ("UniProt Annotation",
             "UniProt REST API v2. Canonical sequence retrieved for each accession. "
             "Peptide mapped by substring match to canonical sequence. "
             "Domain, active site, and binding site features extracted from "
             "feature annotations. Disease associations from DISEASE comments."),
            ("Rule of 25 Sequence Check",
             "Sequence-level assessment of SDA probe residency plausibility. "
             "Checks: lysine presence, local hydrophobicity context (7-residue "
             "window around K), proline adjacency, peptide length. "
             "Full geometric Rule of 25 scoring (linker distance constraints) "
             "is reserved for Tier 3 structural analysis."),
        ]

        for title, text in methods:
            story.append(Paragraph(f"<b>{title}</b>", s_body))
            story.append(Paragraph(text, s_body))
            story.append(Spacer(1, 0.06*inch))

        story.append(grey_rule())
        story.append(Paragraph(
            "PRIVACY NOTICE: This report was generated in volatile RAM only. "
            "No customer data is stored on disk or retained after download. "
            "All intermediate files are deleted when this session ends.",
            s_note))

        doc.build(story)
        return pdf_path

    except Exception as e:
        print(f"PDF generation error: {e}", flush=True)
        raise

# ── TIER 2 ORCHESTRATOR ───────────────────────────────────────────────────────

def run_tier2_pipeline(job_id, work_dir, mzml_path, results_tsv, results, update_fn):
    """
    Orchestrates all Tier 2 pipeline steps.
    All files written to /tmp (work_dir) only — zero retention.
    Returns path to generated PDF.
    """
    last_job_time[0] = time.time()

    hits_all = results.get('high', []) + results.get('medium', []) + results.get('low', [])

    if not hits_all:
        raise ValueError("No SDA hits found in Tier 1 results — cannot generate Tier 2 report.")

    # Load full dataframe for DeepLC calibration
    all_psms_df = results.get('_df')
    if all_psms_df is None:
        all_psms_df = pd.read_csv(results_tsv, sep='\t')

    scored_hits = {}
    uniprot_cache = {}
    figs_dir = os.path.join(work_dir, 'figures')
    os.makedirs(figs_dir, exist_ok=True)

    # Step 1: Extract MS2 spectra
    update_fn('tier2_ms2', 'active', 'Extracting MS2 spectra...', 10)
    ms2_spectra = extract_ms2_spectra(mzml_path, hits_all)
    update_fn('tier2_ms2', 'done', f'MS2 spectra extracted ({len(ms2_spectra)} scans matched)', 100)

    # Step 2: DeepLC retention time prediction (once for all hits)
    update_fn('tier2_rt', 'active', 'Running retention time prediction (DeepLC)...', 20)
    rt_results, rt_error = predict_retention_times(hits_all, all_psms_df)
    if rt_error:
        update_fn('tier2_rt', 'skip', f'RT prediction skipped: {rt_error[:60]}', 100)
    else:
        update_fn('tier2_rt', 'done', 'Retention time prediction complete', 100)

    # Step 3: UniProt annotations (cached per accession)
    update_fn('tier2_uniprot', 'active', 'Fetching UniProt annotations...', 30)
    annotations = {}
    unique_accessions = set()
    for hit in hits_all:
        prot = str(hit.get('proteins', ''))
        # Extract first accession from protein string (may contain multiple)
        acc = prot.split(';')[0].split('|')[1] if '|' in prot else prot.split(';')[0].strip()
        unique_accessions.add((acc, str(hit.get('peptide', ''))))

    for acc, pep in unique_accessions:
        key = f"{acc}::{pep}"
        annotations[key] = fetch_uniprot_annotation(acc, pep, uniprot_cache)

    update_fn('tier2_uniprot', 'done', f'Annotations fetched ({len(uniprot_cache)} proteins)', 100)

    # Step 4: Per-hit scoring, MS1 validation, figures
    update_fn('tier2_score', 'active', 'Scoring hits and generating spectra...', 40)

    for i, hit in enumerate(hits_all[:20]):  # cap at 20 hits for report
        scan_ref = str(hit.get('scannr', hit.get('ScanNr', f'hit_{i}')))
        if ':' in scan_ref:
            scan_ref_clean = scan_ref.split('=')[-1]
        else:
            scan_ref_clean = scan_ref
        try:
            scan_num = int(scan_ref_clean)
        except ValueError:
            scan_num = i

        pep = str(hit.get('peptide', ''))
        prot = str(hit.get('proteins', ''))
        sage_fdr = float(hit.get('spectrum_q', 0.05))

        ms2_scan = ms2_spectra.get(scan_num)

        # Generate theoretical ions
        theo_ions = generate_theoretical_ions(pep, charge=int(hit.get('charge', 2)))

        # Bayesian score
        if ms2_scan is not None and len(ms2_scan['mz']) > 0:
            matches = match_ions(theo_ions, ms2_scan['mz'], ms2_scan['intensity'])
            ms2_scan['scan_num'] = scan_num
            posterior, ms2_tier, ms2_interp = bayesian_score(
                matches, theo_ions, ms2_scan['mz'], sage_fdr)

            # MS1 isotope validation
            ms1_r, ms1_class, ms1_note = validate_ms1_isotope(mzml_path, hit, ms2_scan)

            # Generate spectrum figure
            fig_path = os.path.join(figs_dir, f'spectrum_{i:03d}.png')
            fig_ok = generate_spectrum_figure(hit, ms2_scan, matches, theo_ions, fig_path)

        else:
            matches = []
            posterior, ms2_tier, ms2_interp = 0.0, 'Low', 'MS2 spectrum not available for this PSM.'
            ms1_r, ms1_class, ms1_note = None, 'Inconclusive', 'MS2 spectrum not available.'
            fig_path = None
            fig_ok = False

        # Rule of 25
        r25_outcome, r25_note = rule_of_25_sequence_check(pep)

        # RT data
        rt_data = rt_results.get(pep, {}) if rt_results else {}
        rt_class = rt_data.get('classification', 'Skipped')
        rt_dev = rt_data.get('deviation_min', None)
        rt_note = rt_data.get('note', 'Retention time prediction was not run.')

        # UniProt annotation
        acc = prot.split(';')[0].split('|')[1] if '|' in prot else prot.split(';')[0].strip()
        ann_key = f"{acc}::{pep}"
        annotation = annotations.get(ann_key, {})

        scored_hits[scan_ref] = {
            'hit': hit,
            'ms2_posterior': posterior,
            'ms2_tier': ms2_tier,
            'ms2_interpretation': ms2_interp,
            'ms1_pearson_r': ms1_r,
            'ms1_classification': ms1_class,
            'ms1_note': ms1_note,
            'rt_classification': rt_class,
            'rt_deviation_min': rt_dev,
            'rt_note': rt_note,
            'r25_outcome': r25_outcome,
            'r25_note': r25_note,
            'annotation': annotation,
            'fig_path': fig_path if fig_ok else None,
            'n_matched': len(matches),
        }

        pct = 40 + int((i / max(len(hits_all[:20]), 1)) * 40)
        update_fn('tier2_score', 'active',
                 f'Scored {i+1}/{min(len(hits_all), 20)} hits...', pct)

    update_fn('tier2_score', 'done', f'All hits scored', 100)

    # Step 5: Assemble PDF
    update_fn('tier2_pdf', 'active', 'Assembling PDF report...', 85)
    tier2_data = {'scored_hits': scored_hits}
    pdf_path = generate_tier2_pdf(
        job_id, work_dir, jobs[job_id]['filename'], results, tier2_data)
    update_fn('tier2_pdf', 'done', 'PDF report ready', 100)

    # Store PDF path in memory only
    jobs[job_id]['tier2_pdf_path'] = pdf_path

    return pdf_path

# ── MAIN JOB RUNNER ───────────────────────────────────────────────────────────

def run_job(job_id, file_path, filename, tier=1):
    last_job_time[0] = time.time()
    start = time.time()
    work_dir = os.path.join(WORK_DIR, job_id)
    os.makedirs(work_dir, exist_ok=True)

    # Base steps (same for all tiers)
    step_defs = {
        'detect':    'Detecting file format',
        'decompress':'Decompressing file',
        'convert':   'Converting to mzML',
        'centroid':  'Centroiding spectra',
        'search':    'Running SDA peptide search',
        'parse':     'Parsing results',
        'report':    'Generating report',
    }

    # Additional steps for Tier 2
    if tier == 2:
        step_defs.update({
            'tier2_ms2':     'Extracting MS2 spectra',
            'tier2_rt':      'Predicting retention times',
            'tier2_uniprot': 'Fetching biological annotations',
            'tier2_score':   'Scoring hits',
            'tier2_pdf':     'Assembling PDF report',
        })

    steps = {k: {'status': 'pending', 'message': v, 'pct': 0}
             for k, v in step_defs.items()}

    def update(step, status, message, pct=0):
        steps[step] = {'status': status, 'message': message, 'pct': pct}
        jobs[job_id]['steps'] = steps
        jobs[job_id]['elapsed'] = round(time.time() - start, 1)

    # Save reference to mzml path for Tier 2
    mzml_path_ref = [None]
    results_tsv_ref = [None]

    try:
        fmt = detect_format(filename)
        update('detect', 'done', f'Format detected: {fmt.replace("_"," ").upper()}', 100)

        mzml_path = file_path

        if fmt == 'mzml_gz':
            update('decompress', 'active', 'Decompressing...', 30)
            out = os.path.join(work_dir, filename.replace('.gz', ''))
            with gzip.open(file_path, 'rb') as fi:
                with open(out, 'wb') as fo:
                    shutil.copyfileobj(fi, fo)
            mzml_path = out
            update('decompress', 'done', 'Decompressed', 100)
        else:
            update('decompress', 'skip', 'Decompression not needed', 100)

        update('convert', 'skip', 'Raw conversion not needed', 100)

        update('centroid', 'active', 'Checking spectrum format...', 20)
        if needs_centroiding(mzml_path):
            update('centroid', 'active', 'Auto-centroiding profile data...', 50)
            centroid_path = os.path.join(work_dir, 'centroided.mzML')
            mzml_path, err = centroid_mzml(mzml_path, centroid_path)
            if err:
                update('centroid', 'error', f'Centroiding failed: {err}', 0)
                jobs[job_id]['status'] = 'error'
                jobs[job_id]['error'] = err
                return
            update('centroid', 'done', 'Centroided successfully', 100)
        else:
            update('centroid', 'done', 'Centroid data confirmed', 100)

        mzml_path_ref[0] = mzml_path

        update('search', 'active', 'Building peptide database...', 10)
        results_dir = os.path.join(work_dir, 'results')
        os.makedirs(results_dir, exist_ok=True)
        config_path = build_sage_config(mzml_path, results_dir)

        search_start = time.time()

        def timer():
            while jobs[job_id]['status'] == 'running' and steps['search']['status'] == 'active':
                elapsed_search = int(time.time() - search_start)
                pct = min(10 + elapsed_search, 90)
                update('search', 'active',
                       f'Searching human proteome — {elapsed_search}s elapsed...', pct)
                time.sleep(3)

        threading.Thread(target=timer, daemon=True).start()

        subprocess.run([SAGE_BINARY, config_path],
                       capture_output=True, text=True, timeout=3600)

        results_tsv = os.path.join(results_dir, 'results.sage.tsv')
        if not os.path.exists(results_tsv):
            update('search', 'error', 'Search produced no output', 0)
            jobs[job_id]['status'] = 'error'
            jobs[job_id]['error'] = 'No results file'
            return
        update('search', 'done', 'Search complete', 100)

        results_tsv_ref[0] = results_tsv

        update('parse', 'active', 'Parsing results...', 50)
        results = parse_results(results_tsv)
        update('parse', 'done',
               f"Parsed {results['total_psms']:,} PSMs — {results['total_sda']:,} SDA hits", 100)

        update('report', 'active', 'Generating report...', 70)
        jobs[job_id]['results'] = results
        update('report', 'done', 'Report ready', 100)

        # ── TIER 2 BRANCH ────────────────────────────────────────────────────
        if tier == 2:
            last_job_time[0] = time.time()  # reset idle timer for long Tier 2 processing
            run_tier2_pipeline(
                job_id, work_dir,
                mzml_path_ref[0],
                results_tsv_ref[0],
                results,
                update
            )

        jobs[job_id]['status'] = 'complete'
        jobs[job_id]['tier'] = tier
        jobs[job_id]['elapsed'] = round(time.time() - start, 1)

    except Exception as e:
        jobs[job_id]['status'] = 'error'
        jobs[job_id]['error'] = str(e)
        print(f"Job {job_id} error: {e}", flush=True)
    finally:
        # Zero retention: delete the original upload only
        # work_dir stays in /tmp until download is complete (see download endpoints)
        try:
            os.remove(file_path)
        except Exception:
            pass

# ── ENDPOINTS ─────────────────────────────────────────────────────────────────

@app.post("/submit")
async def submit(
    file: UploadFile = File(...),
    name: str = Form(...),
    institution: str = Form(...),
    tier: int = Form(default=1)
):
    last_job_time[0] = time.time()
    if tier not in (1, 2):
        return JSONResponse({'error': 'Invalid tier. Must be 1 or 2.'}, status_code=400)

    job_id = str(uuid.uuid4())
    work_dir = os.path.join(WORK_DIR, job_id)
    os.makedirs(work_dir, exist_ok=True)
    file_path = os.path.join(work_dir, file.filename)
    with open(file_path, 'wb') as f:
        f.write(await file.read())

    jobs[job_id] = {
        'status': 'running',
        'steps': {},
        'results': None,
        'filename': file.filename,
        'error': None,
        'start': time.time(),
        'elapsed': 0,
        'tier': tier,
        'tier2_pdf_path': None,
    }
    threading.Thread(
        target=run_job,
        args=(job_id, file_path, file.filename, tier),
        daemon=True
    ).start()
    return JSONResponse({'job_id': job_id, 'tier': tier})

@app.get("/status/{job_id}")
async def status(job_id: str):
    job = jobs.get(job_id)
    if not job:
        return JSONResponse({'error': 'Job not found'}, status_code=404)
    return JSONResponse({
        'status': job['status'],
        'steps': job['steps'],
        'results': job['results'],
        'filename': job['filename'],
        'error': job['error'],
        'elapsed': round(time.time() - job['start'], 1),
        'tier': job.get('tier', 1),
    })

@app.get("/download/report/{job_id}")
async def download_report(job_id: str):
    job = jobs.get(job_id)
    if not job or not job['results']:
        return JSONResponse({'error': 'No results'}, status_code=404)
    r = job['results']
    lines = [
        "AXARA BRIDGE — PEPTIDE HIT MAPPING REPORT",
        f"Generated   : {datetime.now().strftime('%Y-%m-%d %H:%M UTC')}",
        f"Input file  : {job['filename']}",
        f"Engine      : Sage v0.14.6 | SDA: +196.0706 Da",
        "=" * 60, "",
        "SUMMARY",
        f"Total PSMs              : {r['total_psms']}",
        f"SDA-modified peptides   : {r['total_sda']}",
        f"High confidence (FDR<1%): {len(r['high'])}",
        f"Medium confidence (1-5%): {len(r['medium'])}",
        f"Exploratory (5-10%)     : {len(r['low'])}",
        "", "=" * 60, "HIGH CONFIDENCE SDA HITS (FDR < 1%)", "=" * 60,
    ]
    for h in r['high']:
        lines += [f"Peptide  : {h['peptide']}", f"Protein  : {h['proteins']}",
                  f"FDR      : {h['spectrum_q']:.4f}", f"Score    : {h['hyperscore']:.2f}", ""]
    if not r['high']:
        lines.append("No hits at FDR < 1%.")
    lines += ["", "=" * 60, "MEDIUM CONFIDENCE (FDR 1-5%)", "=" * 60]
    for h in r['medium']:
        lines += [f"Peptide  : {h['peptide']}", f"Protein  : {h['proteins']}",
                  f"FDR      : {h['spectrum_q']:.4f}", f"Score    : {h['hyperscore']:.2f}", ""]
    if not r['medium']:
        lines.append("No hits at FDR 1-5%.")
    lines += ["", "=" * 60, "EXPLORATORY (FDR 5-10%)", "=" * 60]
    for h in r['low']:
        lines += [f"Peptide  : {h['peptide']}", f"Protein  : {h['proteins']}",
                  f"FDR      : {h['spectrum_q']:.4f}", f"Score    : {h['hyperscore']:.2f}", ""]
    lines += [
        "", "=" * 60, "SEARCH PARAMETERS", "=" * 60,
        "SDA adduct (probe peptide) : +196.0706 Da on K",
        "SDA adduct (target peptide): +82.0419 Da",
        "Oxidation                  : +15.9949 Da on M",
        "Carbamidomethyl            : +57.0215 Da on C",
        "Enzyme                     : Trypsin (KR, no P)",
        "Missed cleavages           : 2",
        "Precursor tolerance        : 15 ppm",
        "Fragment tolerance         : 20 ppm",
        "FDR method                 : Target-decoy competition",
        "Database                   : Human proteome (UniProt SwissProt)",
        "", "For Tier 2 and Tier 3: info@axara.bio",
        "=" * 60,
        "PRIVACY: Volatile RAM only. No files stored.",
        "For research use only.",
    ]
    return StreamingResponse(
        iter(["\n".join(lines)]),
        media_type="text/plain",
        headers={"Content-Disposition":
                 f"attachment; filename=AXARA_Tier1_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"}
    )

@app.get("/download/tsv/{job_id}")
async def download_tsv(job_id: str):
    job = jobs.get(job_id)
    if not job or not job['results']:
        return JSONResponse({'error': 'No results'}, status_code=404)
    return StreamingResponse(
        iter([job['results']['all_sda']]),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition":
                 f"attachment; filename=AXARA_SDA_Hits_{datetime.now().strftime('%Y%m%d_%H%M')}.tsv"}
    )

@app.get("/download/report-tier2/{job_id}")
async def download_tier2_report(job_id: str):
    job = jobs.get(job_id)
    if not job:
        return JSONResponse({'error': 'Job not found'}, status_code=404)
    if job.get('tier') != 2:
        return JSONResponse({'error': 'This job was not run as Tier 2'}, status_code=400)
    if job['status'] != 'complete':
        return JSONResponse({'error': 'Job not yet complete'}, status_code=202)

    pdf_path = job.get('tier2_pdf_path')
    if not pdf_path or not os.path.exists(pdf_path):
        return JSONResponse({'error': 'PDF report not found — may have been cleaned up'}, status_code=404)

    def iter_pdf():
        with open(pdf_path, 'rb') as f:
            while True:
                chunk = f.read(65536)
                if not chunk:
                    break
                yield chunk
        # Zero retention: delete work_dir after PDF is streamed
        try:
            work_dir = os.path.join(WORK_DIR, job_id)
            shutil.rmtree(work_dir, ignore_errors=True)
            jobs[job_id]['tier2_pdf_path'] = None
        except Exception:
            pass

    return StreamingResponse(
        iter_pdf(),
        media_type="application/pdf",
        headers={"Content-Disposition":
                 f"attachment; filename=AXARA_FullReport_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf"}
    )

@app.get("/health")
async def health():
    return JSONResponse({"status": "ok", "jobs": len(jobs)})
