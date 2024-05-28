from distutils.util import execute

from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt


def encode_patient_data(qc, patient_data, mutation_positions):
    binary_data = ""
    for i, nucleotide in enumerate(patient_data):
        if nucleotide == 'A':
            binary_data += '00'
        elif nucleotide == 'C':
            binary_data += '01'
        elif nucleotide == 'G':
            binary_data += '10'
        elif nucleotide == 'T':
            binary_data += '11'

        # Проверка за мутации в гените BRCA1 и BRCA2
        if i in mutation_positions:
            qc.z(i * 2)  # Прилага се Z гейт, ако е открита мутация

    for i, bit in enumerate(binary_data):
        if bit == '1':
            qc.x(i)  # Прилага се X (NOT) гейт за битовете, които трябва да бъдат в състояние |1⟩

    for i in range(len(binary_data)):
        qc.h(i)  # Прилага се Hadamard гейт на всеки кубит за създаване на суперпозиция

    # Допълнителни гейтове
    for i in range(len(binary_data)):
        qc.y(i)  # Добавяне на Y гейт
        qc.z(i)  # Добавяне на Z гейт

    for control_qubit in range(len(binary_data)):
        for target_qubit in range(control_qubit + 1, len(binary_data)):
            qc.cp(0.25 * 3.14159, control_qubit, target_qubit)  # Контролиран фазов гейт


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
    # Примерни данни за пациенти с рак и здрави пациенти
    cancer_patient_data = "ATCGATCG"
    healthy_patient_data = "GCTAGCTA"

    # Позиции на мутациите в гените BRCA1 и BRCA2
    brca1_positions = [2, 5]  # Примерни позиции
    brca2_positions = [3, 6]  # Примерни позиции

    # Уверете се, че данните на двата пациента имат еднаква дължина
    if len(cancer_patient_data) != len(healthy_patient_data):
        raise ValueError("Данните за пациентите с рак и здравите пациенти трябва да са с еднаква дължина")

    # Уверете се, че дължината на данните на пациента е четна
    if len(cancer_patient_data) % 2 != 0:
        raise ValueError("Дължината на данните на пациента трябва да е четна")

    # Определяне на броя кубити
    n_qubits = len(cancer_patient_data) * 2  # Всеки нуклеотид се представя с 2 бита

    # Създаване на квантови вериги за раковите и здравите клетки
    qc_cancer_cells = QuantumCircuit(n_qubits)
    qc_healthy_cells = QuantumCircuit(n_qubits)

    # Кодиране на данните на пациента в квантови състояния
    encode_patient_data(qc_cancer_cells, cancer_patient_data, brca1_positions + brca2_positions)
    encode_patient_data(qc_healthy_cells, healthy_patient_data, brca1_positions + brca2_positions)

    # Прилагане на CNOT гейтове
    for control_qubit in range(n_qubits):
        for target_qubit in range(control_qubit + 1, n_qubits):
            qc_cancer_cells.cx(control_qubit, target_qubit)
            qc_healthy_cells.cx(control_qubit, target_qubit)

    # Измерване на резултатите
    qc_cancer_cells.measure_all()
    qc_healthy_cells.measure_all()

    # Визуализация на квантовите вериги с измервания
    print("Квантова верига за рак:")
    print(qc_cancer_cells.draw())

    print("Квантова верига за здраве:")
    print(qc_healthy_cells.draw())

    # Симулация
    simulator = Aer.get_backend('aer_simulator')

    # Изпълнение на квантовата верига за пациент с рак
    compiled_cancer_circuit = transpile(qc_cancer_cells, simulator)
    result_cancer = execute(compiled_cancer_circuit, simulator).result()
    counts_cancer = result_cancer.get_counts(compiled_cancer_circuit)
    normalized_counts_cancer = {state: count / sum(counts_cancer.values()) for state, count in counts_cancer.items()}

    # Изпълнение на квантовата верига за здрав пациент
    compiled_healthy_circuit = transpile(qc_healthy_cells, simulator)
    result_healthy = execute(compiled_healthy_circuit, simulator).result()
    counts_healthy = result_healthy.get_counts(compiled_healthy_circuit)
    normalized_counts_healthy = {state: count / sum(counts_healthy.values()) for state, count in counts_healthy.items()}

    # Печат на нормализираните резултати за дебъгване
    print("Нормализирани резултати (Рак):", normalized_counts_cancer)
    print("Нормализирани резултати (Здраве):", normalized_counts_healthy)

    # Изчисляване на вероятностите
    cancer_probability, healthy_probability = calculate_probability(normalized_counts_cancer, normalized_counts_healthy)

    # Печат на изчислените вероятности за дебъгване
    print("Вероятност за рак:", cancer_probability)
    print("Вероятност за здраве:", healthy_probability)

    # Сравнение на вероятностите и вземане на решение
    result = compare_probabilities(cancer_probability, healthy_probability)
    print("Резултат от сравнението:", result)

    # Визуализация на резултатите
    plot_histogram([counts_cancer, counts_healthy], legend=["Cancer", "Healthy"])
    plt.show()


if __name__ == "__main__":
    main()