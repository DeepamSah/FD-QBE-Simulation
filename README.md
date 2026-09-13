# FD-QBE Simulation

Repository for the simulation and experimental validation of the Forward-Deniable Quantum Ballot Encryption (FD-QBE) protocol.

## Contents

### Source Code
- fdqbe_simulation.py

### Experimental Data
- results.csv
- security_metrics.csv
- monte_carlo_results.csv
- shamir_results.csv
- delay_safety_results.csv

### Figures
- Encryption time analysis
- Throughput analysis
- Memory usage analysis
- Delay safety analysis

## Main Findings

- Near-zero mutual information between votes and ciphertexts.
- Information-theoretic forward deniability.
- Delay-independent ballot secrecy.
- Efficient ballot encryption and tallying.

## Installation

```bash
pip install -r requirements.txt
```

## Run

```bash
python fdqbe_simulation.py
```

## Citation

If you use this repository in academic work, please cite:

Sah, D., Dev, B., and Kshetri, K.
FD-QBE Simulation Repository.
GitHub Repository.
2026.

## Author

Deepam Sah,
Bibhan Dev,
Kshitiz Kshetri*

## Contact
Email: kkshetri16@gmail.com*, 
      sahdeepam12@gmail.com
