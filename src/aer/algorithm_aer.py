from qiskit import QuantumCircuit, Aer, transpile, assemble, execute
from qiskit.providers.aer import Aer

from qiskit.visualization import plot_histogram

import matplotlib.pyplot as plt


def encode_patient_data(qc, patient_data):
    binary_data = ""

    for nucleotide in patient_data:
        if nucleotide == 'A':
            binary_data += '00'
        elif nucleotide == 'C':
            binary_data += '01'
        elif nucleotide == 'G':
            binary_data += '10'
        elif nucleotide == 'T':
            binary_data += '11'

    for i, bit in enumerate(binary_data):
        if bit == '1':
            qc.x(i)  # Apply X (NOT) gate for bits that should be in state |1⟩

    for i in range(len(binary_data)):
        qc.h(i)  # Apply Hadamard gate to each qubit to create superposition


def calculate_probability(normalized_counts_cancer, normalized_counts_healthy):
    total_cancer = sum(normalized_counts_cancer.values())
    total_healthy = sum(normalized_counts_healthy.values())
    total_samples = total_cancer + total_healthy

    cancer_probability = (total_cancer / total_samples) * 100
    healthy_probability = (total_healthy / total_samples) * 100

    return cancer_probability, healthy_probability


def compare_probabilities(cancer_probability, healthy_probability, threshold=50):
    if cancer_probability > threshold and cancer_probability > healthy_probability:
        return "Cancer"

    elif healthy_probability > threshold and healthy_probability > cancer_probability:
        return "Healthy"

    else:
        return "Undetermined"


def main():
    # Example data for cancer and healthy patients
    cancer_patient_data = "ATCGATCG"
    healthy_patient_data = "GCTAGCTA"

    # Ensure both patient data have the same length

    if len(cancer_patient_data) != len(healthy_patient_data):
        raise ValueError("Cancer and healthy patient data must be of the same length")

    # Creating quantum circuits for cancer and healthy cells
    n_qubits = len(cancer_patient_data) * 2  # Each nucleotide is represented by 2 bits
    qc_cancer_cells = QuantumCircuit(n_qubits)
    qc_healthy_cells = QuantumCircuit(n_qubits)

    # Encoding patient data into quantum states
    encode_patient_data(qc_cancer_cells, cancer_patient_data)
    encode_patient_data(qc_healthy_cells, healthy_patient_data)

    # Apply CX (CNOT) gates
    for control_qubit in range(n_qubits):
        for target_qubit in range(control_qubit + 1, n_qubits):
            qc_cancer_cells.cx(control_qubit, target_qubit)
            qc_healthy_cells.cx(control_qubit, target_qubit)

    # Measurement of results
    qc_cancer_cells.measure_all()
    qc_healthy_cells.measure_all()

    # Visualization of quantum circuits with measurements
    print("Quantum cancer circuit:")
    print(qc_cancer_cells.draw())
    print("Quantum healthy circuit:")
    print(qc_healthy_cells.draw())

    # Simulation
    simulator = Aer.get_backend('aer_simulator')

    # Execute cancer patient circuit
    compiled_cancer_circuit = transpile(qc_cancer_cells, simulator)
    result_cancer = execute(compiled_cancer_circuit, simulator).result()
    counts_cancer = result_cancer.get_counts(qc_cancer_cells)
    normalized_counts_cancer = {state: count / sum(counts_cancer.values()) for state, count in counts_cancer.items()}

    # Execute healthy patient circuit
    compiled_healthy_circuit = transpile(qc_healthy_cells, simulator)
    result_healthy = execute(compiled_healthy_circuit, simulator).result()
    counts_healthy = result_healthy.get_counts(qc_healthy_cells)
    normalized_counts_healthy = {state: count / sum(counts_healthy.values()) for state, count in counts_healthy.items()}

    # Calculate probabilities
    cancer_probability, healthy_probability = calculate_probability(normalized_counts_cancer, normalized_counts_healthy)

    # Compare probabilities and make a decision
    result = compare_probabilities(cancer_probability, healthy_probability)

    print("Comparison result:", result)

    # Visualization of results
    plot_histogram([counts_cancer, counts_healthy], legend=["Cancer", "Healthy"])
    plt.show()


if __name__ == "__main__":
    main()