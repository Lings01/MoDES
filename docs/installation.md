# Installation

## Requirements

- Python >= 3.10
- pip

## Quick Install

```bash
git clone https://github.com/Lings01/MoDES.git
cd MoDES
pip install -e .
```

## Development Install

```bash
pip install -e .
pip install -r requirements-dev.txt
```

## Verify

```bash
python -c "from modes import MoDES, MoDEData; print('OK')"
python -m pytest -q
```

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| numpy | >=1.21 | Numerical arrays |
| scipy | >=1.7 | Statistical distributions |
| pandas | >=1.3 | Data frames |
| statsmodels | >=0.13 | NB GLM fitting |
| anndata | >=0.8 | Single-cell data |
| matplotlib | >=3.5 | Plotting |
| seaborn | >=0.11 | Heatmaps |
| networkx | >=2.6 | Graph export |

### Optional

| Package | Version | Purpose |
|---|---|---|
| pytest | >=7.0 | Testing |
| pytest-cov | >=4.0 | Coverage |
