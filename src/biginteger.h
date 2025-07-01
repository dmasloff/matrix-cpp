//
// Created by Даниил Маслов on 29.11.2022.
//
#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class BigInteger;

class Rational;

BigInteger operator+(BigInteger a, const BigInteger& b);

BigInteger operator-(BigInteger a, const BigInteger& b);

BigInteger operator*(BigInteger a, const BigInteger& b);

BigInteger operator/(BigInteger a, const BigInteger& b);

BigInteger operator%(BigInteger a, const BigInteger& b);

bool operator==(const BigInteger& a, const BigInteger& b);

bool operator!=(const BigInteger& a, const BigInteger& b);

bool operator<(const BigInteger& a, const BigInteger& b);

bool operator>(const BigInteger& a, const BigInteger& b);

bool operator<=(const BigInteger& a, const BigInteger& b);

bool operator>=(const BigInteger& a, const BigInteger& b);

std::ostream& operator<<(std::ostream& os, const BigInteger& bigInteger);

std::istream& operator>>(std::istream& is, BigInteger& bigInteger);

BigInteger operator""_bi(unsigned long long num);

BigInteger operator""_bi(const char* c, size_t len);

Rational operator+(Rational a, const Rational& b);

Rational operator-(Rational a, const Rational& b);

Rational operator*(Rational a, const Rational& b);

Rational operator/(Rational a, const Rational& b);

bool operator==(const Rational& a, const Rational& b);

bool operator!=(const Rational& a, const Rational& b);

bool operator<(const Rational& a, const Rational& b);

bool operator>(const Rational& a, const Rational& b);

bool operator<=(const Rational& a, const Rational& b);

bool operator>=(const Rational& a, const Rational& b);

std::ostream& operator<<(std::ostream& os, const Rational& r);

class BigInteger {
  public:
    BigInteger(long long num);

    BigInteger(const std::string& str);

    BigInteger();

    BigInteger& operator=(BigInteger other);

    BigInteger operator+() const;

    BigInteger operator-() const;

    BigInteger& operator+=(const BigInteger& other);

    BigInteger& operator-=(const BigInteger& other);

    BigInteger& operator*=(long long num);

    BigInteger& operator*=(const BigInteger& other);

    BigInteger& operator/=(const BigInteger& other);

    BigInteger& operator%=(const BigInteger& other);

    BigInteger& operator++();

    BigInteger operator++(int);

    BigInteger& operator--();

    BigInteger operator--(int);

    std::string toString() const;

    bool isNegative() const;

    void changeSign();

    void makePositive();

    void Swap(BigInteger& other);

    explicit operator bool() const;

    explicit operator int() const;

    ~BigInteger();

  private:
    static inline const long long kBase_ = 1e9;
    std::vector<long long> digits_;
    long long size_;
    bool isNegative_;

    long long operator[](long long ind) const;

    void shrinkDigits();

    void subtractionWithConstSign(const BigInteger& other);

    void subtractionWithChangeableSign(const BigInteger& other);

    void addition(const BigInteger& other);

    void commonCaseAddition(const BigInteger& other, bool isAddition);

    std::pair<BigInteger, BigInteger> division(const BigInteger& other);

    bool absGreater(const BigInteger& b) const;

    friend bool operator==(const BigInteger& a, const BigInteger& b);

    friend bool operator<(const BigInteger& a, const BigInteger& b);
};

BigInteger gcd(const BigInteger& a, const BigInteger& b);

class Rational {
  public:
    Rational(int num);

    Rational(const BigInteger& bi);

    Rational();

    Rational& operator=(const Rational& r);

    Rational& operator+=(const Rational& r);

    Rational& operator-=(const Rational& r);

    Rational& operator*=(const Rational& r);

    Rational& operator/=(const Rational& r);

    Rational operator-() const;

    Rational operator+() const;

    explicit operator double() const;

    std::string toString() const;

    std::string asDecimal(size_t precision) const;

    ~Rational();

    friend bool operator==(const Rational& a, const Rational& b);

    friend bool operator<(const Rational& a, const Rational& b);

  private:
    template <typename Field>
    friend struct is_field;
    BigInteger numerator_, denominator_;

    void normalize();
};

//
// Residue
//

template <size_t L, size_t R, size_t N>
struct is_greater {
    static const bool value = ((L + R) / 2) * ((L + R) / 2) > N;
};

template <size_t L, size_t R, size_t N>
static const bool is_greater_v = is_greater<L, R, N>::value;

template <size_t L, size_t R, size_t N>
struct square_root {
    static const size_t value =
        (R == L + 1 ? L
                    : square_root<(is_greater_v<L, R, N> ? L : (L + R) / 2),
                                  (is_greater_v<L, R, N> ? (L + R) / 2 : R),
                                  N>::value);
};

template <size_t N>
const static size_t square_root_v = square_root<0, N, N>::value;

template <size_t I, size_t N>
struct is_prime {
    static const bool value = (N % I != 0 && is_prime<I - 1, N>::value);
};

template <size_t N>
struct is_prime<1, N> {
    static const bool value = true;
};

template <size_t N>
static const bool is_prime_v = is_prime<square_root_v<N>, N>::value;

template <>
static const bool is_prime_v<0> = false;

template <>
static const bool is_prime_v<1> = false;

template <typename Field>
struct is_field {
    static const bool isField = Field::kIsPrime_;
};

template <>
struct is_field<Rational> {
    static const bool isField = true;
};

template <typename Field>
static const bool is_field_v = is_field<Field>::isField;

template <size_t N>
class Residue {
  public:
    static const bool kIsPrime_ = is_prime_v<N>;

    Residue(const Residue& another) : value(another.value){};

    Residue(int value) : value(value < 0 ? N - ((-value) % N) : value % N){};

    Residue() : Residue(0) {}

    Residue operator+=(const Residue& another) {
        value += another.value;
        value %= N;
        return *this;
    }

    Residue operator-=(const Residue& another) {
        value += (value < another.value ? N : 0);
        value -= another.value;
        return *this;
    }

    Residue operator*=(const Residue& another) {
        value *= another.value;
        value %= N;
        return *this;
    }

    Residue operator/=(const Residue& another) {
        static_assert(kIsPrime_, "Division in not a prime field!");
        *this *= another.inverseElement();
        return *this;
    }

    Residue operator=(const Residue& another) {
        value = another.value;
        return value;
    }

    explicit operator int() const {
        return value;
    }

    bool operator==(const Residue& another) const {
        return value == another.value;
    }

    size_t characteristic() const {
        return N;
    }

  private:
    size_t value = 0;

    Residue inverseElement() const {
        size_t power = N - 2;
        if (power == 0) {
            return Residue(1);
        }
        Residue ans = 1, multiply = *this;
        while (power > 0) {
            if (power % 2 == 1) {
                ans *= multiply;
                --power;
            } else {
                multiply *= multiply;
                power /= 2;
            }
        }
        return ans;
    }
};

template <size_t N>
Residue<N> operator+(Residue<N> a, const Residue<N>& b) {
    return a += b;
}

template <size_t N>
Residue<N> operator-(Residue<N> a, const Residue<N>& b) {
    return a -= b;
}

template <size_t N>
Residue<N> operator*(Residue<N> a, const Residue<N>& b) {
    return a *= b;
}

template <size_t N>
Residue<N> operator/(Residue<N> a, const Residue<N>& b) {
    return a /= b;
}

template <size_t N>
std::ostream& operator<<(std::ostream& os, const Residue<N>& another) {
    os << static_cast<int>(another);
    return os;
}

#endif