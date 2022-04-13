#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
class Matrix {
public:
    int m_rows, m_cols;
    double *m_values;
    Matrix(); // По умолчанию создается пустая матрица
    Matrix(int rows, int cols); // Выделяет память, элементы не инициализирует
    Matrix(const Matrix& other); // Копирование матрицы
    ~Matrix(); //Деструктор

    bool isEqual(const Matrix& other);
    Matrix matrixPow(int n); // Возведение матрицы в степень
    Matrix matrixCpy(); // Копирование матрицы
    bool isZero(); // Нулевая ли
    bool isUnit(); // Еденичная ли
    bool isDiagonal(); // Диагональная ли
    bool isSymmetrical(); // Симметрическая ли
    bool isUpperTriangular(); // Верхнетреугольная ли
    bool isLowerTriangular(); // Нижнетреугольная ли
    void checkType(); // Проверка типа матрицы

    bool isMultipliable(const Matrix& other); // Умножаемые ли
    double get(int i, int j) const; // Возврат элемента по индексу
    bool readFromKeyboard(); // Чтение матрицы с клавиатуры
    void printToScreen(); // Печать в консоль
    int rows() const; // Возврат числа строк матрицы
    int cols() const; // Возврат числа столбцов матрицы
    Matrix transpose(); // Транспонирование матрицы
    void resize(int rows, int cols); // Переопределение размеров матрицы
    void multiply(double value); // Умножение матрицы на число

    int choicefunction(int k, int K); // Вспомогательная для cutmatrix функция
    Matrix cutmatrix(int I, int J); // Алгебраическое дополнение матрицы
    double elemDet(); // Определитель матрицы 2x2
    double Det(); // Определитель произвольной квадратной матрицы
    Matrix inversematrix(); // Обратная матрица

    Matrix operator+=(const Matrix& other);
    Matrix operator-=(const Matrix& other);
    Matrix operator*=(const Matrix& other);
};
//Создание класса-наследника "Vector"
class Vector : public Matrix {
public:
    int length;
    Vector();
    Vector(int len);
    Vector(const Vector& other);

    double norm(); //Норма вектора
    void printVector(); //Печать в консоль
    double scalarProduct(const Vector& other); //Скалярное произведение векторов
    Vector vectorProduct(const Vector& other); //Векторное произведение векторов
};
Vector::Vector() : Matrix() {
    length = 0;
}
Vector::Vector(int len) : Matrix(len, 1) {
    length = len;
}
Vector::Vector(const Vector& other) : Matrix(other){
    length = other.length;
}

double Vector::norm() {
    double result = 0;
    for(int i = 0; i < length; i++) {
        result += m_values[i] * m_values[i];
    }
    return result;
}
double Vector::scalarProduct(const Vector& other) {
    double result = 0;
    for(int i = 0; i < length; i++) {
        result += m_values[i] * other.m_values[i];
    }
    return result;
}
Vector Vector::vectorProduct(const Vector& other) {
    double x, y, z;
    if(length != 3 || other.length != 3) {
        cout << "Error: Vector product is not defined for that type of vectors" << endl;
    }
    x = m_values[1] * other.m_values[2] - other.m_values[1] * m_values[2];
    y = other.m_values[0] * m_values[2] - m_values[0] * other.m_values[2];
    z = m_values[0] * other.m_values[1] - other.m_values[0] * m_values[1];
    m_values[0] = x;
    m_values[1] = y;
    m_values[2] = z;
    return *this;
}
void Vector::printVector() {
    cout << "(";
    for (int i = 0; i < length - 1; i++) {
        cout << m_values[i] << ", ";
    }
    cout << m_values[length - 1] << ")" << endl;
}
//--------------------------------------------------------------------------------------------
Matrix::Matrix() {
    m_rows = 0;
    m_rows = 0;
    m_values = nullptr;
}
Matrix::Matrix(int rows, int cols) {
    m_rows = rows;
    m_cols = cols;
    m_values = new double[m_cols * m_rows];
}
Matrix::Matrix(const Matrix& other) {
    m_values = new double[other.m_cols * other.m_rows];
    for(int i = 0; i < other.m_rows * other.m_cols; i++) {
        m_values[i] = other.m_values[i];
    }
}
Matrix::~Matrix() {
    if (nullptr != m_values)
        delete[] m_values;
}

bool Matrix::isEqual(const Matrix& other) {
    if(m_rows != other.m_rows || m_cols != other.m_cols) {
        return false;
    }
    for(int i = 0; i < m_rows * m_cols; i++) {
        if(m_values[i] != other.m_values[i]) return false;
    }
    return true;
}
Matrix Matrix::matrixCpy() {
    Matrix result(m_rows, m_cols);
    for(int i = 0; i < m_rows * m_cols; i++) {
        result.m_values[i] = m_values[i];
    }
    return result;
}
Matrix Matrix::matrixPow(int n) {
    Matrix result = (*this).matrixCpy();
    for(int i = 0; i < n - 1; i++) {
        result *= (*this);
    }
    for(int i = 0; i < m_rows * m_cols; i++) {
        m_values[i] = result.m_values[i];
    }
    return *this;
}
bool Matrix::isZero() {
    for(int i = 0; i < m_rows * m_cols; i++) {
        if(m_values[i] != 0) return false;
    }
    return true;
}
bool Matrix::isUnit() {
    if(m_rows != m_cols) return false;
    for(int i = 0; i < m_rows * m_cols; i++) {
        if(i % (m_rows + 1) == 0) {
            if(m_values[i] != 1) return false;
        }
        else {
            if(m_values[i] != 0) return false;
        }
    }
    return true;
}
bool Matrix::isDiagonal() {
    if(m_rows != m_cols) return false;
    for(int i = 0; i < m_rows * m_cols; i++) {
        if((i % (m_rows + 1) != 0) && (m_values[i] != 0)) return false;
    }
    return true;
}
bool Matrix::isSymmetrical() {
    bool boolean;
    if(m_rows != m_cols) return false;
    Matrix tmp = (*this).matrixCpy();
    tmp.transpose();
    boolean = (*this).isEqual(tmp);
    tmp.~Matrix();
    return boolean;
}
bool Matrix::isUpperTriangular() {
    if(m_rows != m_cols) return false;
    for(int i = 1; i < m_rows; i++) {
        for(int j = 0; j < i; j++) {
            if(m_values[i*m_cols + j] != 0) return false;
        }
    }
    return true;
}
bool Matrix::isLowerTriangular() {
    int boolean;
    if(m_rows != m_cols) return false;
    Matrix tmp = (*this).matrixCpy();
    tmp.transpose();
    boolean = tmp.isUpperTriangular();
    tmp.~Matrix();
    return boolean;
}
void Matrix::checkType() {
    if(m_rows != m_cols) {
        cout << "Matrix Type: <Non-square matrix>" << endl;
        return void();
    }
    if((*this).isZero()) {
        cout << "Matrix Type: <Zero matrix>" << endl;
        return void();
    }
    if((*this).isUnit()) {
        cout << "Matrix Type: <Unit matrix>" << endl;
        return void();
    }
    if((*this).isDiagonal()) {
        cout << "Matrix Type: <Diagonal matrix>" << endl;
        return void();
    }
    if((*this).isUpperTriangular()) {
        cout << "Matrix Type: <UpperTriangular matrix>" << endl;
        return void();
    }
    if((*this).isLowerTriangular()) {
        cout << "Matrix Type: <LowerTriangular matrix>" << endl;
        return void();
    }
    cout << "Matrix Type: <Square matrix>" << endl;
}

double Matrix::get(int i, int j) const {
    return m_values[i*m_cols + j];
}
bool Matrix::isMultipliable(const Matrix& other) {
    return m_cols == other.m_rows;
}
bool Matrix::readFromKeyboard() {
    if(m_values == nullptr) {
        cout << "Error: No memory has been allocated";
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < m_rows; i++) {
        for(int j = 0; j < m_cols; j++) {
            cin >> m_values[i*m_cols + j];
        }
    }
    return true;
}
void Matrix::printToScreen() {
    for(int i = 0; i < m_cols; i++) cout << "------";
    cout << endl;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            cout << fixed << setprecision(2) << m_values[i*m_cols + j] << " ";
        }
        cout << endl;
    }
    for(int i = 0; i < m_cols; i++) cout << "------";
    cout << endl;
}
Matrix Matrix::transpose() {
    Matrix other(m_cols, m_rows);
    for(int i = 0; i < m_rows; i++) {
        for(int j = 0; j < m_cols; j++) {
            other.m_values[j*other.m_cols + i] = m_values[i*m_cols + j];
        }
    }
    this->resize(m_cols, m_rows);
    for(int i = 0; i < m_rows * m_cols; i++){
        m_values[i] = other.m_values[i];
    }
    return *this;
}
int Matrix::rows() const {
    return m_rows;
}
int Matrix::cols() const {
    return m_cols;
}
void Matrix::resize(int rows, int cols) {
    m_rows = rows;
    m_cols = cols;
    if (m_values)
        delete[] m_values;
    m_values = new double[m_rows * m_cols];
}
void Matrix::multiply(double value) {
    for(int i = 0; i < m_rows * m_cols; i++) {
        m_values[i] *= value;
    }
}
int Matrix::choicefunction(int k, int K) {
    return k - K > 0 ? 1 : 0;
}
Matrix Matrix::cutmatrix(int I, int J) {
    Matrix result(m_rows - 1, m_rows - 1);
    int iflag = 0;
    int jflag = 0;
    for(int i = 0; i < m_rows; i++) {
        for(int j = 0; j < m_rows; j++) {
            if(j == J) jflag = 1;
            if(i == I) iflag = 1;
            if(j != J && i != I) {
                result.m_values[(i - iflag*choicefunction(i, I)) * (m_rows - 1) + j - jflag*choicefunction(j, J)] = m_values[i * m_rows + j];
            }
        }
    }
    return result;
}
double Matrix::elemDet() {
    if(m_rows != m_cols)
    {
        cout << "Error: Determinant is not defined for that type of matrix" << endl;
        exit(EXIT_FAILURE);
    }
    if(m_rows != 2)
    {
        cout << "Error: Calculating elemDet for not 2x2 matrix" << endl;
        exit(EXIT_FAILURE);
    }
    return (m_values[0] * m_values[3] - m_values[1] * m_values[2]);
}
double Matrix::Det() {
    if(m_rows != m_cols) {
        cout << "Error: Determinant is not defined for that type of matrix" << endl;
        exit(EXIT_FAILURE);
    }
    if(m_rows == 1) {
        return m_values[0];
    }
    if(m_rows == 2) {
        return (*this).elemDet();
    }
    double det = 0;
    for(int i = 0; i < m_rows; i++) {
        det += pow((-1), i) * (m_values[i] * (*this).cutmatrix(0, i).Det());
    }
    return det;
}
Matrix Matrix::inversematrix() {
    if(m_rows != m_cols) {
        cout << "Error: Determinant is not defined for that type of matrix" << endl;
        exit(EXIT_FAILURE);
    }
    double det;
    det = (*this).Det();
    if(det == 0) {
        cout << "Error: Zero determinant. The inverse matrix doesn`t exist" << endl;
        exit(EXIT_FAILURE);
    }
    Matrix result(m_rows, m_cols);
    for(int i = 0; i < m_rows; i++) {
        for(int j = 0; j < m_cols; j++) {
           result.m_values[i*m_cols + j] = pow((-1), i+j) * (*this).cutmatrix(i, j).Det();
        }
    }
    result.transpose();
    result.multiply(1 / det);
    this->resize(m_cols, m_rows);
    for(int i = 0; i < m_rows * m_cols; i++) {
        m_values[i] = result.m_values[i];
    }
    result.~Matrix();
    return *this;
}

Matrix Matrix::operator+= (const Matrix& other) {
    if(m_rows != other.m_rows || m_cols != other.m_cols){
        cout << "Error: Addition of matrices of different dimensions" << endl;
        exit(EXIT_FAILURE);
    }
    for(size_t i = 0; i < m_rows * m_cols; i++) {
        m_values[i] = m_values[i] + other.m_values[i];
    }
    return *this;
}
Matrix Matrix::operator-= (const Matrix& other) {
    if(m_rows != other.m_rows || m_cols != other.m_cols){
        cout << "Error: Addition of matrices of different dimensions" << endl;
        exit(EXIT_FAILURE);
    }
    for(size_t i = 0; i < m_rows * m_cols; i++) {
        m_values[i] = m_values[i] - other.m_values[i];
    }
    return *this;
}
Matrix Matrix::operator*= (const Matrix& other){
    if (!isMultipliable(other)) {
        cout << "Error: Matrices aren`t isMultipliable";
        exit(EXIT_FAILURE);
    }
    Matrix result(m_rows, other.m_cols);
    for(int i = 0; i < m_rows; i++) {
        for(int j = 0; j < other.m_cols; j++) {
            result.m_values[i*other.m_cols + j] = 0;
            for(int k = 0; k < m_cols; k++){
                result.m_values[i*other.m_cols + j] += get(i, k) * other.get(k, j);
            }
        }
    }
    this->resize(m_rows, other.m_cols);
    for(int i = 0; i < m_rows * other.m_cols; i++) {
        m_values[i] = result.m_values[i];
    }
    result.~Matrix();
    return *this;
}
int main() {
    int work_mode; //Режим работы программы (Матричный или Векторный)
    int switcher; // Переключатель

    cout << "Welcome to the matrix calculator!" << endl;
    cout << "Please select the operating mode of the program:" << endl;
    cout << "1)Matrix mode" << endl;
    cout << "2)Vector mode" << endl;
    cin >> work_mode;
    switch(work_mode){
        case(1):
        {
            int default_rows; // Номер строк первоначальной матрицы
            int default_cols; // Номер столбцов первоначальной матрицы
            double mul_num; // Переменная для case(4), число, на которое умножается матрица
            int exp_num; //Переменная для case(9), степень, в которую возводим матрицу

            cout << "Enter the size of the matrix:" << endl;
            cout << "Rows: ";
            cin >> default_rows;
            cout << "Columns: ";
            cin >> default_cols;
            Matrix M(default_rows, default_cols);
            cout << "Please enter the matrix elements line by line:" << endl;
            M.readFromKeyboard();
            /*
            cout << M.isZero() << " Zero?" << endl;
            cout << M.isUnit() << " Unit?" << endl;
            cout << M.isDiagonal() << " Diag?" << endl;
            cout << M.isSymmetrical() << " Symm?" << endl;
            cout << M.isUpperTriangular() << " UT?" << endl;
            cout << M.isLowerTriangular() << " LT?" << endl; 
            */
            cout << "Input completed successfully." << endl;
            cout << "Please, select the operation you want to perform:" << endl;
            while(1) {
                cout << endl;
                cout << "0)Change the initial matrix" << endl;
                cout << "1)Add a matrix" << endl;
                cout << "2)Subtract the matrix" << endl;
                cout << "3)Multiply by the matrix" << endl;
                cout << "4)Multiply by a number" << endl;
                cout << "5)Transpose the matrix" << endl;
                cout << "6)Inverse the matrix" << endl;
                cout << "7)Determinant of the matrix" << endl;
                cout << "8)Print current matrix" << endl;
                cout << "9)Exponentiation of the matrix" << endl;
                cout << "10)Check matrix type" << endl;
                cout << "11)Exit" << endl;


                cin >> switcher;
                switch(switcher) {
                    case(0):
                    {
                        cout << "Enter the size of the matrix:" << endl;
                        cout << "Rows: ";
                        cin >> default_rows;
                        cout << "Columns: ";
                        cin >> default_cols;
                        M.resize(default_rows, default_cols);
                        cout << "Please enter the matrix elements line by line:" << endl;
                        M.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Please, select the operation you want to perform:" << endl;
                        break;
                    }
                    case(1):
                    {
                        cout << "Enter the matrix to be added line by line:" << endl;
                        cout << "Remember that the size of your matrix should be (" << M.rows() << ", " << M.cols() << ")" << endl;
                        Matrix m(default_rows, default_cols);
                        m.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Result:" << endl;
                        M += m;
                        M.printToScreen();
                        m.~Matrix();
                        break;
                    }
                    case(2):
                    {
                        cout << "Enter the matrix to be subtracted line by line:" << endl;
                        cout << "Remember that the size of your matrix should be (" << M.rows() << ", " << M.cols() << ")" << endl;
                        Matrix m(default_rows, default_cols);
                        m.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Result:" << endl;
                        M -= m;
                        M.printToScreen();
                        m.~Matrix();
                        break;
                    }
                    case(3):
                    {
                        cout << "Enter the size of the matrix:" << endl;
                        cout << "Remember that the size of your matrix should be (" << M.cols() << ", " << "any" << ")" << endl;
                        cout << "Rows: ";
                        cin >> default_rows;
                        cout << "Columns: ";
                        cin >> default_cols;
                        cout << "Enter the matrix to be multiplied line by line:" << endl;
                        Matrix m(default_rows, default_cols);
                        m.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Result:" << endl;
                        M *= m;
                        M.printToScreen();
                        m.~Matrix();
                        break;
                    }
                    case(4):
                    {
                        cout << "Enter the number to be multiplied:" << endl;
                        cin >> mul_num;
                        cout << "Input completed successfully." << endl;
                        cout << "Result:" << endl;
                        M.multiply(mul_num);
                        M.printToScreen();
                        break;                
                    }
                    case(5):
                    {
                        M.transpose();
                        M.printToScreen();
                        break;
                    }
                    case(6):
                    {
                        if(M.rows() != M.cols()) {
                            cout << "Non-square matrix. The inverse matrix is not defined" << endl;
                            break;
                        }
                        if(M.Det() == 0) {
                            cout << "Zero determinant. The inverse matrix doesn`t exist" << endl;
                            break;
                        }
                        M.inversematrix();
                        M.printToScreen();
                        break;
                    }
                    case(7):
                    {
                        if(M.rows() != M.cols()) {
                            cout << "Non-square matrix. The determinant is not defined" << endl;
                            break;
                        }
                        cout << "Matrix Determinant = " << M.Det() << endl;
                        break;
                    }
                    case(8):
                    {
                        M.printToScreen();
                        break;
                    }
                    case(9):
                    {
                            cout << "Enter the power:" << endl;
                            cout << "Remember that the number must be natural" << endl;
                            cin >> exp_num;
                            cout << "Input completed successfully." << endl;
                            cout << "Result:" << endl;
                            M.matrixPow(exp_num);
                            M.printToScreen();
                            break;
                    }
                    case(10):
                    {
                        M.checkType();
                        break;
                    }
                    case(11):
                    {
                        exit(EXIT_SUCCESS);
                    }
                    default:
                    {
                        cout << "Error: Wrong input" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        case(2):
        {
            int default_length;
            cout << "Enter the length of the vector:" << endl;
            cout << "Length: ";
            cin >> default_length;
            Vector V(default_length);
            cout << "Please enter the vector elements:" << endl;
            V.readFromKeyboard();
            cout << "Input completed successfully." << endl;
            cout << "Please, select the operation you want to perform:" << endl;
            while(1) {
                cout << endl;
                cout << "0)Change the initial vector" << endl;
                cout << "1)Scalar product" << endl;
                cout << "2)Vector product" << endl;
                cout << "3)Print current vector" << endl;
                cout << "4)Exit" << endl;

                cin >> switcher;
                switch(switcher) {
                    case(0):
                    {
                        cout << "Enter the length of the vector:" << endl;
                        cout << "Length: ";
                        cin >> default_length;
                        V.resize(default_length, 1);
                        cout << "Please enter the vector elements:" << endl;
                        V.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Please, select the operation you want to perform:" << endl;
                        break;
                    }
                    case(1):
                    {
                        cout << "Enter the vector with which you want to calculate the scalar product:" << endl;
                        Vector v(default_length);
                        v.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Result: " << V.scalarProduct(v) << endl;
                        v.~Matrix();
                        break;
                    }
                    case(2):
                    {
                        if(default_length != 3) {
                            cout << "Vector product is not defined for a vector of this dimension, it is defined only for 3d vectors" << endl;
                            break;
                        }
                        cout << "Enter the vector with which you want to calculate the vector product:" << endl;
                        Vector v(default_length);
                        v.readFromKeyboard();
                        cout << "Input completed successfully." << endl;
                        cout << "Result:" << endl;
                        V.vectorProduct(v);
                        V.printVector();
                        v.~Matrix();
                        break;
                    }
                    case(3):
                    {
                        V.printVector();
                        break;
                    }
                    case(4):
                    {
                        exit(EXIT_SUCCESS);
                    }
                    default:
                    {
                        cout << "Error: Wrong input" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        default:
        {
            cout << "Error: Wrong input" << endl;
            exit(EXIT_FAILURE);   
        }
    }
    
}