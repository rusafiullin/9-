import numpy as np
import matplotlib.pyplot as plt
from math import pow, log
import sys


class Nuclear_Analyzer:
    def __init__(self):
        self.a1 = 15.7
        self.a2 = 17.3
        self.a3 = 0.71
        self.a4 = 23.6
        self.a5 = 34

        self.mp = 1.007
        self.mn = 1.008
        self.me = 0.0005

    def _validate_input(self, Z, A):
        if not isinstance(Z, (int, float)) or not isinstance(A, (int, float)):
            raise ValueError("Z и A должны быть числами")

        Z = int(Z)
        A = int(A)

        if Z < 0:
            raise ValueError("Z (атомный номер) не может быть отрицательным")
        if A < 0:
            raise ValueError("A (массовое число) не может быть отрицательным")
        if Z > A:
            raise ValueError("Z не может быть больше A")
        if A == 0 and Z != 0:
            raise ValueError("При A=0 Z также должен быть 0")

        return Z, A

    def formula(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            raise e

        if A == 0:
            return 0

        try:
            first = self.a1 * A
            second = -self.a2 * (A ** (2 / 3))
            third = -self.a3 * (Z ** 2 / (A ** (1 / 3)))
            four = -self.a4 * ((A - 2 * Z) ** 2 / A)

            if Z % 2 == 0 and (A - Z) % 2 == 0:
                five = self.a5 * (A ** (-3 / 4))
            elif Z % 2 == 1 and (A - Z) % 2 == 1:
                five = -self.a5 * (A ** (-3 / 4))
            else:
                five = 0

            B = first + second + third + four + five
            return B

        except (ZeroDivisionError, OverflowError) as e:
            raise ValueError(f"Ошибка вычисления для Z={Z}, A={A}: {e}")

    def ud_energy(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            raise e

        if A == 0:
            return 0

        B = self.formula(Z, A)
        if B == 0:
            return 0

        return B / A

    def mass(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            raise e

        if A == 0:
            return 0

        try:
            B = self.formula(Z, A)
            M = Z * self.mp + (A - Z) * self.mn - B / 931.49
            return max(0, M)  # Масса не может быть отрицательной

        except (ValueError, TypeError) as e:
            raise ValueError(f"Ошибка вычисления массы для Z={Z}, A={A}: {e}")

    def radius(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            raise e

        if A == 0:
            return 0

        try:
            r0 = 1.2
            radius = r0 * (A ** (1 / 3))
            return max(0, radius)  # Радиус не может быть отрицательным

        except (ValueError, OverflowError) as e:
            raise ValueError(f"Ошибка вычисления радиуса для Z={Z}, A={A}: {e}")

    def beta(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            return f"Ошибка: {e}"

        if A == 0:
            return "Невозможно определить"

        try:
            mc = self.mass(Z, A)

            if Z + 1 <= A:
                mplus = self.mass(Z + 1, A)
            else:
                mplus = float('inf')

            if Z - 1 >= 0:
                mminus = self.mass(Z - 1, A)
            else:
                mminus = float('inf')

            betaminus = mc > mplus and mplus != float('inf')
            betaplus = mc > mminus and mminus != float('inf')

            if not betaminus and not betaplus:
                return "Устойчивое"
            elif betaminus:
                return "Неустойчивое бета-минус"
            else:
                return "Неустойчивое бета-плюс"

        except Exception as e:
            return f"Ошибка определения бета-распада: {e}"

    def delenie(self, Z, A):
        try:
            Z, A = self._validate_input(Z, A)
        except ValueError as e:
            return f"Ошибка: {e}"

        if A == 0:
            return "Невозможно определить"

        try:
            if A % 2 != 0 or Z % 2 != 0:
                return "Нечетное A или Z"

            A_osk = A // 2
            Z_osk = Z // 2

            if A_osk % 2 != 0 or Z_osk % 2 != 0:
                return "Осколки нечетные"
            else:
                return "Возможно"

        except Exception as e:
            return f"Ошибка проверки деления: {e}"


def safe_input(prompt, input_type=float, min_val=None, max_val=None):
    while True:
        try:
            value = input_type(input(prompt))
            if min_val is not None and value < min_val:
                print(f"Значение не может быть меньше {min_val}")
                continue
            if max_val is not None and value > max_val:
                print(f"Значение не может быть больше {max_val}")
                continue
            return value
        except ValueError:
            print("Пожалуйста, введите корректное число")
        except KeyboardInterrupt:
            print("\nПрограмма прервана пользователем")
            sys.exit(0)


def interactive_mode():
    analyzer = Nuclear_Analyzer()

    print("\n=== ИНТЕРАКТИВНЫЙ РЕЖИМ АНАЛИЗА ЯДЕР ===")

    while True:
        print("\n" + "=" * 50)
        print("Введите параметры ядра (или 'q' для выхода):")

        try:
            Z_input = input("Атомный номер Z: ").strip()
            if Z_input.lower() == 'q':
                break

            A_input = input("Массовое число A: ").strip()
            if A_input.lower() == 'q':
                break

            Z = int(Z_input)
            A = int(A_input)
            specific_energy = analyzer.ud_energy(Z, A)
            mass = analyzer.mass(Z, A)
            radius = analyzer.radius(Z, A)
            beta_stable = analyzer.beta(Z, A)
            fission = analyzer.delenie(Z, A)

            print(f"\nРезультаты анализа ядра Z={Z}, A={A}:")
            print(f"  Удельная энергия связи: {specific_energy:.2f} МэВ/нуклон")
            print(f"  Масса атома: {mass:.4f} а.е.м.")
            print(f"  Радиус ядра: {radius:.2f} фм")
            print(f"  Бета-распад: {beta_stable}")
            print(f"  Деление: {fission}")

        except ValueError as e:
            print(f"Ошибка ввода: {e}")
        except Exception as e:
            print(f"Неожиданная ошибка: {e}")


def graphics():
    analyzer = Nuclear_Analyzer()

    isotopes = [
        ("U-238", 92, 238),
        ("Pu-239", 94, 239),
        ("Cf-252", 98, 252),
        ("Pu-238", 94, 238),
        ("Te-135", 52, 135),
        ("Ni-60", 28, 60),
        ("O-16", 8, 16),
        ("N-15", 7, 15),
        ("P-29", 15, 29),
        ("Si-29", 14, 29),
        ("Cr-52", 24, 52)
    ]

    print("АНАЛИЗ АТОМНЫХ ЯДЕР")
    results = []

    for name, Z, A in isotopes:
        try:
            specific_energy = analyzer.ud_energy(Z, A)
            mass = analyzer.mass(Z, A)
            radius = analyzer.radius(Z, A)
            beta_stable = analyzer.beta(Z, A)
            fission = analyzer.delenie(Z, A)

            results.append({
                'name': name,
                'Z': Z,
                'A': A,
                'specific_energy': specific_energy,
                'mass': mass,
                'radius': radius,
                'beta_stable': beta_stable,
                'fission': fission
            })

            print(f"\n{name} (Z={Z}, A={A}):")
            print(f"  Удельная энергия связи: {specific_energy:.2f} МэВ/нуклон")
            print(f"  Масса атома: {mass:.4f} а.е.м.")
            print(f"  Радиус ядра: {radius:.2f} фм")
            print(f"  Бета-распад: {beta_stable}")
            print(f"  Деление: {fission}")

        except Exception as e:
            print(f"Ошибка анализа {name}: {e}")
            continue

    return results, isotopes


def create_plots(results, isotopes):
    if not results:
        print("Нет данных для построения графиков")
        return

    try:
        plt.figure(figsize=(15, 5))

        # График 1: Радиус vs A
        plt.subplot(1, 3, 1)
        A_values = [res['A'] for res in results]
        radii = [res['radius'] for res in results]
        names = [res['name'] for res in results]

        plt.scatter(A_values, radii, color='blue', s=50)
        for i, name in enumerate(names):
            plt.annotate(name, (A_values[i], radii[i]), xytext=(5, 5),
                         textcoords='offset points', fontsize=8)

        A_range = np.linspace(min(A_values), max(A_values), 100)
        r0 = 1.2
        radius_theory = r0 * np.power(A_range, 1 / 3)
        plt.plot(A_range, radius_theory, 'r--', alpha=0.7, label='R = 1.2 × A¹ᐟ³')

        plt.xlabel('Массовое число A')
        plt.ylabel('Радиус ядра (фм)')
        plt.title('Радиус атомного ядра')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.subplot(1, 3, 2)
        specific_energies = [res['specific_energy'] for res in results]

        plt.scatter(A_values, specific_energies, color='green', s=50)
        for i, name in enumerate(names):
            plt.annotate(name, (A_values[i], specific_energies[i]), xytext=(5, 5),
                         textcoords='offset points', fontsize=8)

        plt.xlabel('Массовое число A')
        plt.ylabel('Удельная энергия связи (МэВ/нуклон)')
        plt.title('Удельная энергия связи')
        plt.grid(True, alpha=0.3)

        plt.subplot(1, 3, 3)
        Z_values = [res['Z'] for res in results]

        plt.scatter(Z_values, specific_energies, color='red', s=50)
        for i, name in enumerate(names):
            plt.annotate(name, (Z_values[i], specific_energies[i]), xytext=(5, 5),
                         textcoords='offset points', fontsize=8)

        plt.xlabel('Атомный номер Z')
        plt.ylabel('Удельная энергия связи (МэВ/нуклон)')
        plt.title('Удельная энергия связи vs Z')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Ошибка построения графиков: {e}")


def main():
    print("АНАЛИЗАТОР АТОМНЫХ ЯДЕР")
    print("1 - Интерактивный режим")
    print("2 - Анализ предопределенных изотопов")
    print("3 - Выход")

    while True:
        choice = input("\nВыберите режим работы (1-3): ").strip()

        if choice == '1':
            interactive_mode()
        elif choice == '2':
            try:
                results, isotopes = graphics()
                create_plots(results, isotopes)
            except Exception as e:
                print(f"Ошибка в режиме анализа: {e}")
        elif choice == '3':
            print("Выход из программы")
            break
        else:
            print("Неверный выбор. Пожалуйста, введите 1, 2 или 3")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nПрограмма прервана пользователем")
    except Exception as e:
        print(f"Критическая ошибка: {e}")