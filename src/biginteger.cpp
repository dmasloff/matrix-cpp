#include "biginteger.h"

BigInteger::BigInteger(long long num)
    : digits_(0), size_(0), isNegative_(num < 0) {
    num = (isNegative_ ? -num : num);
    while (num != 0) {
        digits_.push_back(num % kBase_);
        num /= kBase_;
        ++size_;
    }
}

BigInteger::BigInteger(const std::string &str)
    : digits_(0), size_(0), isNegative_(str[0] == '-') {
    long long i = str.size();
    while (i >= 9) {
        digits_.push_back(stoll(str.substr(i - 9, 9)));
        i -= 9;
        ++size_;
    }
    if (isNegative_ && i > 1) {
        digits_.push_back(stoll(str.substr(1, i - 1)));
        ++size_;
    }
    if (!isNegative_ && i > 0) {
        digits_.push_back(stoll(str.substr(0, i)));
        ++size_;
    }
}

BigInteger::BigInteger() : digits_(0), size_(0), isNegative_(false) {}

BigInteger &BigInteger::operator=(BigInteger other) {
    Swap(other);
    return *this;
}

BigInteger BigInteger::operator+() const {
    BigInteger res = *this;
    return res;
}

BigInteger BigInteger::operator-() const {
    BigInteger res = *this;
    if (res != 0) {
        res.isNegative_ = !res.isNegative_;
    }
    return res;
}

BigInteger &BigInteger::operator+=(const BigInteger &other) {
    commonCaseAddition(other, true);
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &other) {
    commonCaseAddition(other, false);
    return *this;
}

BigInteger &BigInteger::operator*=(long long num) {
    isNegative_ = isNegative_ != (num < 0);
    num = (num < 0 ? -num : num);
    long long carry = 0;
    for (long long i = 0; i < size_; ++i) {
        digits_[i] *= num;
        digits_[i] += carry;
        carry = digits_[i] / kBase_;
        digits_[i] %= kBase_;
    }
    while (carry != 0) {
        digits_.push_back(carry % kBase_);
        carry /= kBase_;
        ++size_;
    }
    shrinkDigits();
    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &other) {
    bool signcopy = isNegative_ != other.isNegative_;
    BigInteger copy = 0, tmp;
    isNegative_ = false;
    std::vector<long long> zero;
    for (long long i = 0; i < other.size_; ++i) {
        tmp = (*this);
        tmp *= other.digits_[i];
        tmp.digits_.insert(tmp.digits_.begin(), zero.begin(), zero.end());
        tmp.size_ += i;
        zero.push_back(0);
        copy += tmp;
    }
    *this = copy;
    isNegative_ = signcopy;
    shrinkDigits();
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &other) {
    *this = division(other).first;
    shrinkDigits();
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &other) {
    *this = division(other).second;
    shrinkDigits();
    return *this;
}

BigInteger &BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger res = *this;
    *this += 1;
    return res;
}

BigInteger &BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator--(int) {
    BigInteger res = *this;
    *this -= 1;
    return res;
}

std::string BigInteger::toString() const {
    std::string str, tmp_str;
    if (isNegative_) {
        str.push_back('-');
    }
    long long ind = size_ - 1;
    str += std::to_string((*this)[ind--]);
    while (ind >= 0) {
        tmp_str = std::to_string(digits_[ind]);
        for (size_t i = tmp_str.size(); i < 9; ++i) {
            str.push_back('0');
        }
        str += tmp_str;
        --ind;
    }
    return str;
}

bool BigInteger::isNegative() const {
    return isNegative_;
}

void BigInteger::changeSign() {
    isNegative_ = !isNegative_;
}

void BigInteger::makePositive() {
    isNegative_ = false;
}

void BigInteger::Swap(BigInteger &other) {
    std::swap(this->digits_, other.digits_);
    std::swap(this->size_, other.size_);
    std::swap(this->isNegative_, other.isNegative_);
}

BigInteger::operator bool() const {
    return size_ != 0;
}

BigInteger::operator int() const {
    int res = 0, tmpBase = 1;
    for (long long i = 0; i < size_; ++i) {
        res += digits_[i] * tmpBase;
        tmpBase *= kBase_;
    }
    return res;
}

BigInteger::~BigInteger() = default;

long long BigInteger::operator[](long long ind) const {
    if (0 <= ind && ind < size_) {
        return digits_[ind];
    }
    return 0;
}

void BigInteger::shrinkDigits() {
    while (size_ > 0 && digits_[size_ - 1] == 0) {
        --size_;
    }
    digits_.resize(size_);
    if (size_ == 0) {
        isNegative_ = false;
    }
}

void BigInteger::subtractionWithConstSign(const BigInteger &other) {
    long long carry = 0, i;
    for (i = 0; i < other.size_; ++i) {
        digits_[i] -= other.digits_[i] + carry;
        carry = 0;
        if (digits_[i] < 0) {
            digits_[i] += kBase_;
            carry = 1;
        }
    }
    while (carry > 0) {
        digits_[i] -= carry;
        carry = 0;
        if (digits_[i] < 0) {
            digits_[i] += kBase_;
            carry = 1;
        }
        ++i;
    }
}

void BigInteger::subtractionWithChangeableSign(const BigInteger &other) {
    isNegative_ = !isNegative_;
    size_ = std::max(size_, other.size_);
    digits_.resize(size_, 0);
    long long carry = 0;
    for (long long i = 0; i < other.size_; ++i) {
        digits_[i] = other.digits_[i] - digits_[i] - carry;
        carry = 0;
        if (digits_[i] < 0) {
            digits_[i] += kBase_;
            carry = 1;
        }
    }
}

void BigInteger::addition(const BigInteger &other) {
    long long carry = 0;
    if (size_ < other.size_) {
        digits_.resize(other.size_, 0);
        size_ = other.size_;
    }
    for (long long i = 0; i < size_; ++i) {
        digits_[i] += other[i] + carry;
        carry = digits_[i] / kBase_;
        digits_[i] %= kBase_;
    }
    if (carry > 0) {
        digits_.push_back(carry);
        ++size_;
    }
}

void BigInteger::commonCaseAddition(const BigInteger &other, bool isAddition) {
    if ((isNegative_ != other.isNegative_) == isAddition) {
        if (absGreater(other)) {
            subtractionWithConstSign(other);
        } else {
            subtractionWithChangeableSign(other);
        }
    } else {
        addition(other);
    }
    shrinkDigits();
}

std::pair<BigInteger, BigInteger> BigInteger::division(
    const BigInteger &other) {
    if (size_ < other.size_) {
        return {0, *this};
    }
    bool signcopy = isNegative_ != other.isNegative_;
    long long length = size_ - other.size_;
    isNegative_ = false;
    BigInteger copy, divider = other;
    divider.isNegative_ = false;
    divider.size_ += length;
    divider.digits_.insert(divider.digits_.begin(), size_ - other.size_, 0);
    copy.isNegative_ = false;
    copy.size_ = length + 1;
    copy.digits_.resize(copy.size_, 0);
    for (long long i = length; i >= 0; --i) {
        long long l = 0, r = kBase_, mid;
        BigInteger tmp = divider;
        while (r - l > 1 && size_ != 0) {
            mid = (l + r) / 2;
            tmp = divider;
            tmp *= mid;
            if (!(*this < tmp)) {
                l = mid;
            } else {
                r = mid;
            }
        }
        shrinkDigits();
        copy.digits_[i] = l;
        *this -= divider * l;
        divider.digits_.erase(divider.digits_.begin());
        --divider.size_;
    }
    isNegative_ = signcopy;
    copy.isNegative_ = signcopy;
    shrinkDigits();
    copy.shrinkDigits();
    return {copy, *this};
}

bool BigInteger::absGreater(const BigInteger &b) const {
    if (size_ != b.size_) {
        return size_ > b.size_;
    }
    long long i = size_ - 1;
    while (i >= 0 && digits_[i] == b.digits_[i]) {
        --i;
    }
    return (*this)[i] > b[i];
}

BigInteger operator+(BigInteger a, const BigInteger &b) {
    a += b;
    return a;
}

BigInteger operator-(BigInteger a, const BigInteger &b) {
    a -= b;
    return a;
}

BigInteger operator*(BigInteger a, const BigInteger &b) {
    a *= b;
    return a;
}

BigInteger operator/(BigInteger a, const BigInteger &b) {
    a /= b;
    return a;
}

BigInteger operator%(BigInteger a, const BigInteger &b) {
    a %= b;
    return a;
}

bool operator==(const BigInteger &a, const BigInteger &b) {
    return a.isNegative_ == b.isNegative_ && a.size_ == b.size_ &&
           a.digits_ == b.digits_;
}

bool operator!=(const BigInteger &a, const BigInteger &b) {
    return !(a == b);
}

bool operator<(const BigInteger &a, const BigInteger &b) {
    if (a.isNegative_ != b.isNegative_) {
        return a.isNegative_;
    }
    if (a.size_ != b.size_) {
        return a.isNegative_ != (a.size_ < b.size_);
    }
    long long i = a.size_ - 1;
    while (i >= 0 && a.digits_[i] == b.digits_[i]) {
        --i;
    }
    return a.isNegative_ != (a[i] < b[i]);
}

bool operator>(const BigInteger &a, const BigInteger &b) {
    return (b < a);
}

bool operator<=(const BigInteger &a, const BigInteger &b) {
    return !(a > b);
}

bool operator>=(const BigInteger &a, const BigInteger &b) {
    return !(a < b);
}

std::ostream &operator<<(std::ostream &os, const BigInteger &bigInteger) {
    os << bigInteger.toString();
    return os;
}

std::istream &operator>>(std::istream &is, BigInteger &bigInteger) {
    std::string str;
    is >> str;
    bigInteger = BigInteger(str);
    return is;
}

BigInteger operator ""_bi(unsigned long long num) {
    return BigInteger(num);
}

BigInteger operator ""_bi(const char *c, size_t len) {
    return BigInteger(std::string(c, len));
}

BigInteger gcd(const BigInteger &a, const BigInteger &b) {
    BigInteger acopy = a, bcopy = b;
    acopy.makePositive();
    bcopy.makePositive();
    if (acopy > bcopy) {
        acopy.Swap(bcopy);
    }
    while (acopy != 0) {
        bcopy %= acopy;
        bcopy.Swap(acopy);
    }
    return bcopy;
}

Rational::Rational(int num) : numerator_(num), denominator_(1) {};

Rational::Rational(const BigInteger &bi) : numerator_(bi), denominator_(1) {};

Rational::Rational() : numerator_(0), denominator_(1) {};

Rational &Rational::operator=(const Rational &r) {
    numerator_ = r.numerator_;
    denominator_ = r.denominator_;
    return *this;
}

Rational &Rational::operator+=(const Rational &r) {
    numerator_ *= r.denominator_;
    numerator_ += denominator_ * r.numerator_;
    denominator_ *= r.denominator_;
    normalize();
    return *this;
}

Rational &Rational::operator-=(const Rational &r) {
    numerator_ *= r.denominator_;
    numerator_ -= denominator_ * r.numerator_;
    denominator_ *= r.denominator_;
    normalize();
    return *this;
}

Rational &Rational::operator*=(const Rational &r) {
    numerator_ *= r.numerator_;
    denominator_ *= r.denominator_;
    normalize();
    return *this;
}

Rational &Rational::operator/=(const Rational &r) {
    numerator_ *= r.denominator_;
    denominator_ *= r.numerator_;
    if (denominator_.isNegative()) {
        numerator_.changeSign();
        denominator_.changeSign();
    }
    normalize();
    return *this;
}

Rational Rational::operator-() const {
    Rational res = *this;
    res.numerator_.changeSign();
    return res;
}

Rational Rational::operator+() const {
    return *this;
}

Rational::operator double() const {
    return std::stod(asDecimal(100));
}

std::string Rational::toString() const {
    std::string str;
    str += numerator_.toString();
    if (denominator_ != 1) {
        str += '/';
        str += denominator_.toString();
    }
    return str;
}

std::string Rational::asDecimal(size_t precision) const {
    std::string str;
    BigInteger div = numerator_ / denominator_, mod = numerator_ % denominator_;
    if (div == 0_bi && numerator_.isNegative()) {
        str += '-';
    }
    str += div.toString();
    mod.makePositive();
    if (denominator_ != 1_bi) {
        str.push_back('.');
        mod *= BigInteger("1" + std::string(precision, '0'));
        mod /= denominator_;
        std::string tmp = mod.toString();
        if (tmp.size() < precision) {
            str += std::string(precision - tmp.size(), '0');
        }
        str += tmp;
    }
    return str;
}

Rational::~Rational() = default;

void Rational::normalize() {
    BigInteger tmp = gcd(numerator_, denominator_);
    if (tmp != 1 && tmp != 0) {
        numerator_ /= tmp;
        denominator_ /= tmp;
    };
}

Rational operator+(Rational a, const Rational &b) {
    a += b;
    return a;
}

Rational operator-(Rational a, const Rational &b) {
    a -= b;
    return a;
}

Rational operator*(Rational a, const Rational &b) {
    a *= b;
    return a;
}

Rational operator/(Rational a, const Rational &b) {
    a /= b;
    return a;
}

bool operator==(const Rational &a, const Rational &b) {
    return (a.denominator_ * b.numerator_ == a.numerator_ * b.denominator_);
}

bool operator!=(const Rational &a, const Rational &b) {
    return !(a == b);
}

bool operator<(const Rational &a, const Rational &b) {
    return (a.numerator_ * b.denominator_ < b.numerator_ * a.denominator_);
}

bool operator>(const Rational &a, const Rational &b) {
    return (b < a);
}

bool operator<=(const Rational &a, const Rational &b) {
    return !(a > b);
}

bool operator>=(const Rational &a, const Rational &b) {
    return !(a < b);
}

std::ostream &operator<<(std::ostream &os, const Rational &r) {
    os << r.toString();
    return os;
}

std::istream &operator>>(std::istream &is, Rational &r) {
    int tmp;
    is >> tmp;
    r = Rational(tmp);
    return is;
}