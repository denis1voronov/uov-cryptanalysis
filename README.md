This repository contains Python code to reproduce the experiments of partial key exposure attack (proposed by author) in bit erasure model on UOV schemes.
We use cryptographic-estimators package (https://github.com/Crypto-TII/CryptographicEstimators), licensed under Apache-2.0. To run attack estimator, execute the following:

```bash
git clone https://github.com/denis1voronov/uov-cryptanalysis
cd uov-cryptanalysis
python -m venv .venv

# Windows
.venv\Scripts\activate

# Linux/macOS
source .venv/bin/activate

pip install -r requirements.txt
python estimations.py
```
