from qiskit import QuantumCircuit, transpile, assemble
from qiskit_aer import Aer

from qiskit.visualization import plot_histogram, plot_bloch_multivector

import matplotlib.pyplot as plt


def read_sequence_from_file(file_path):
    with open(file_path, 'r') as file:
        sequence = file.readline().strip()

    return sequence


def align_sequences(seq1, seq2):
    max_len = max(len(seq1), len(seq2))
    padded_seq1 = seq1.ljust(max_len, '-')
    padded_seq2 = seq2.ljust(max_len, '-')

    return padded_seq1, padded_seq2


def encode_patient_data(qc, patient_data, healthy_data):
    binary_data = ""

    mutation_detected = False

    for i, nucleotide in enumerate(patient_data):

        if nucleotide == 'A':

            binary_data += '00'

        elif nucleotide == 'C':

            binary_data += '01'

        elif nucleotide == 'G':

            binary_data += '10'

        elif nucleotide == 'T':

            binary_data += '11'

            # Проверка за мутации чрез сравняване със здравата секвенция

        if i < len(healthy_data) and nucleotide != healthy_data[i]:
            qc.z(i * 2)  # Прилага се Z гейт, ако е открита мутация

            mutation_detected = True

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

    return mutation_detected


def calculate_probability(normalized_counts_cancer, normalized_counts_healthy, mutation_detected):
    total_cancer = sum(normalized_counts_cancer.values())

    total_healthy = sum(normalized_counts_healthy.values())

    total_samples = total_cancer + total_healthy

    cancer_probability = (total_cancer / total_samples) * 100

    healthy_probability = (total_healthy / total_samples) * 100

    # Увеличаване на вероятността за рак, ако е открита мутация

    if mutation_detected:
        cancer_probability += 10  # Примерно увеличение на вероятността с 10%

    return cancer_probability, healthy_probability


def compare_probabilities(cancer_probability, healthy_probability, threshold=50):
    if cancer_probability > threshold and cancer_probability > healthy_probability:

        return "Cancer"

    elif healthy_probability > threshold and healthy_probability > cancer_probability:

        return "Healthy"

    else:

        return "Undetermined"


def main():
    # Въвеждане на пътища до файловете от потребителя
    # cancer_file_path = input("Въведете пътя до файла с данни за пациент с рак: ")
    # healthy_file_path = input("Въведете пътя до файла с данни за здрав пациент: ")
    #
    # # Четене на данни за пациенти с рак и здрави пациенти от файлове
    # cancer_patient_data = read_sequence_from_file(cancer_file_path)
    # healthy_patient_data = read_sequence_from_file(healthy_file_path)

    cancer_patient_data = "ATCGATCG"
    healthy_patient_data = "GCTAGCTA"

    # Подравняване на секвенциите
    cancer_patient_data, healthy_patient_data = align_sequences(cancer_patient_data, healthy_patient_data)

    # Определяне на броя кубити
    n_qubits = len(cancer_patient_data) * 2  # Всеки нуклеотид се представя с 2 бита

    # Създаване на квантови вериги за раковите и здравите клетки
    qc_cancer_cells = QuantumCircuit(n_qubits)
    qc_healthy_cells = QuantumCircuit(n_qubits)

    # Кодиране на данните на пациента в квантови състояния и откриване на мутации
    mutation_detected_cancer = encode_patient_data(qc_cancer_cells, cancer_patient_data, healthy_patient_data)
    mutation_detected_healthy = encode_patient_data(qc_healthy_cells, healthy_patient_data, healthy_patient_data)

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
    # result_cancer_assemble = assemble(compiled_cancer_circuit, simulator)
    result_cancer = simulator.run(compiled_cancer_circuit).result()
    counts_cancer = result_cancer.get_counts(compiled_cancer_circuit)
    normalized_counts_cancer = {state: count / sum(counts_cancer.values()) for state, count in counts_cancer.items()}

    # Изпълнение на квантовата верига за здрав пациент
    compiled_healthy_circuit = transpile(qc_healthy_cells, simulator)
    # result_healthy_assemble = assemble(compiled_healthy_circuit, simulator)
    result_healthy = simulator.run(compiled_healthy_circuit).result()
    counts_healthy = result_healthy.get_counts(compiled_healthy_circuit)
    normalized_counts_healthy = {state: count / sum(counts_healthy.values()) for state, count in counts_healthy.items()}

    # Печат на нормализираните резултати за дебъгване
    print("Нормализирани резултати (Рак):", normalized_counts_cancer)
    print("Нормализирани резултати (Здраве):", normalized_counts_healthy)

    # Изчисляване на вероятностите
    cancer_probability, healthy_probability = calculate_probability(normalized_counts_cancer, normalized_counts_healthy,
                                                                    mutation_detected_cancer)

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
