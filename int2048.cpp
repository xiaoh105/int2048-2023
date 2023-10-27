#include <iostream>
#include <int2048.h>

sjtu::polynomial::polynomial()
{
  len = 1;
  a = new __int128 [1]{0};
}
sjtu::polynomial::polynomial(const sjtu::int2048 &val,
                             bool high_digit_first = false)
{
  len = val.len;
  a = new __int128 [len + 5];
  if (high_digit_first)
    for (int i = 0; i < val.len; ++i) a[i] = val.a[val.len - i - 1];
  else
    for (int i = 0; i < val.len; ++i) a[i] = val.a[i];
}

sjtu::polynomial::polynomial(const sjtu::polynomial &val)
{
  len = val.len;
  a = new __int128 [len + 5];
  for (int i = 0; i < len; ++i) a[i] = val.a[i];
}

sjtu::polynomial::~polynomial()
{
  delete [] a;
}

void sjtu::polynomial::ExtendLen(int new_len)
{
  auto *new_a = new __int128 [new_len + 6];
  for (int i = 0; i < len; ++i) new_a[i] = a[i];
  for (int i = len; i < new_len; ++i) new_a[i] = 0;
  delete [] a;
  a = new_a;
  len = new_len;
}

void sjtu::polynomial::ChangeIndex()
{
  int *rev;
  rev = new int [len + 5];
  rev[0] = 0;
  for (int i = 1; i < len; ++i)
  {
    rev[i] = (rev[i >> 1] >> 1);
    if (i & 1) rev[i] += (len >> 1);
  }
  for (int i = 0; i < len; ++i)
    if (i < rev[i]) std::swap(a[i], a[rev[i]]);
  delete [] rev;
}

__int128 sjtu::pow_mod(__int128 base, __int128 pow)
{
  __int128 ret = 1, cur_pow = base;
  while (pow > 0)
  {
    if (pow & 1)
    {
      ret *= cur_pow;
      ret %= sjtu::polynomial::mod;
    }
    cur_pow *= cur_pow, cur_pow %= sjtu::polynomial::mod;
    pow >>= 1;
  }
  return ret;
}

void sjtu::Extend_GCD(__int128 a, __int128 b, __int128 &x, __int128 &y)
{
  if (b == 0)
  {
    x = 1, y = 0;
    return;
  }
  Extend_GCD(b, a % b, y, x);
  y -= x * (a / b);
}

__int128 sjtu::inverse(__int128 a)
{
  __int128 x, y;
  Extend_GCD(a, sjtu::polynomial::mod, x, y);
  x %= sjtu::polynomial::mod;
  if (x < 0) x += sjtu::polynomial::mod;
  return x;
}

void sjtu::polynomial::NTT(int is_NTT)
{
  ChangeIndex();
  for (int step = 2; step <= len; step <<= 1)
  {
    __int128 w = 0;
    if (is_NTT == 1) { w = pow_mod(root, (mod - 1) / step); }
    else { w = pow_mod(inv, (mod - 1) / step); }
    for (int i = 0; i < len; i += step)
    {
      __int128 cur_w = 1;
      int j = i + (step >> 1);
      for (int k = 0; k < (step >> 1); ++k)
      {
        __int128 f = a[i + k], g = a[j + k];
        a[i + k] = f + cur_w * g % mod, a[i + k] %= mod;
        a[j + k] = f - cur_w * g % mod;
        a[j + k] = (a[j + k] % mod + mod) % mod;
        cur_w *= w, cur_w %= mod;
      }
    }
  }
  if (is_NTT == -1)
  {
    __int128 inv_len = inverse(len);
    for (int i = 0; i < len; ++i) a[i] *= inv_len, a[i] %= mod;
  }
}

sjtu::polynomial &sjtu::polynomial::Multiply(sjtu::polynomial val)
{
  int new_len = 1;
  while (new_len < len * 2 || new_len < val.len * 2) new_len *= 2;
  ExtendLen(new_len);
  val.ExtendLen(new_len);
  NTT(1), val.NTT(1);
  for (int i = 0; i < len; ++i) a[i] *= val.a[i], a[i] %= mod;
  NTT(-1);
  return *this;
}

sjtu::int2048 sjtu::polynomial::ToInteger(bool high_digit_first = false,
                                          int length = 0)
{
  sjtu::int2048 ret;
  delete [] ret.a;
  if (high_digit_first)
  {
    ret.a = new int [length + 5];
    for (int i = 0; i < length; ++i) ret.a[i] = a[length - i - 1];
    while (ret.a[len - 1] == 0 && ret.len >= 2) --ret.len;
    return ret;
  }
  ret.len = len;
  __int128 *tmp;
  tmp = new __int128 [len + 5];
  ret.a = new int [len + 5];
  for (int i = 0; i < len + 5; ++i) ret.a[i] = tmp[i] = 0;
  for (int i = 0; i < len; ++i) tmp[i] = a[i];
  for (int i = 1; i < len; ++i)
  {
    tmp[i] += tmp[i - 1] / sjtu::int2048::base;
    tmp[i - 1] %= sjtu::int2048::base;
  }
  while (tmp[ret.len - 1] >= sjtu::int2048::base)
  {
    tmp[ret.len] += tmp[ret.len - 1] / sjtu::int2048::base;
    tmp[ret.len - 1] %= sjtu::int2048::base;
    ++ret.len;
  }
  while(tmp[ret.len - 1] == 0 && ret.len >= 2) --ret.len;
  for (int i = 0; i < ret.len; ++i) ret.a[i] = tmp[i];
  delete [] tmp;
  return ret;
}

sjtu::int2048::int2048()
{
  len = 1;
  a = new int [1] {0};
  sgn = 1;
}

sjtu::int2048::int2048(long long val)
{
  sgn = 1;
  if (val < 0)
  {
    sgn = -1;
    val = -val;
  }
  len = 0;
  long long tmp = val;
  while (tmp)
  {
    ++len;
    tmp /= base;
  }
  if (val == 0) len = 1; // 0 has length 1
  a = new int [len + 5];
  for (int i = 0; i < len; ++i)
  {
    a[i] = val % base;
    val /= base;
  }
}

sjtu::int2048::int2048(const std::string &s)
{
  len = 0;
  std::string input;
  if (s[0] == '-')
  {
    sgn = -1;
    input = s.substr(1);
    if (input[0] == '0' && input.length() == 1) sgn = 1;
  }
  else
  {
    sgn = 1;
    input = s;
  }
  a = new int [input.length() / base_log10 + 5];
  for (int i = input.length() - 1; i >= 0; i -= base_log10)
  {
    int cur_digit = 0, pow10 = 1;
    for (int j = 0; j < base_log10 && i - j >= 0; ++j)
    {
      cur_digit += pow10 * (input[i - j] - '0');
      pow10 *= 10;
    }
    a[len] = cur_digit, ++len;
  }
}

sjtu::int2048::int2048(const sjtu::int2048 &val)
{
  len = val.len;
  sgn = val.sgn;
  a = new int [len + 5];
  for (int i = 0; i < len; ++i) a[i] = val.a[i];
}

sjtu::int2048::~int2048()
{
  delete [] a;
}

void sjtu::int2048::read(const std::string &s)
{
  len = 0;
  std::string input;
  if (s[0] == '-')
  {
    sgn = -1;
    input = s.substr(1);
    if (input.length() == 1 && input[0] == '0') sgn = 1;
  }
  else
  {
    sgn = 1;
    input = s;
  }
  delete [] a;
  a = new int [input.length() / base_log10 + 5];
  for (int i = input.length() - 1; i >= 0; i -= base_log10)
  {
    int cur_digit = 0, pow10 = 1;
    for (int j = 0; j < base_log10 && i - j >= 0; ++j)
    {
      cur_digit += pow10 * (input[i - j] - '0');
      pow10 *= 10;
    }
    a[len] = cur_digit, ++len;
  }
}

void sjtu::int2048::print()
{
  if (sgn == -1 && (len != 1 || a[0] != 0)) printf("-");
  for (int i = len - 1; i >= 0; --i)
    if (i != len - 1) { printf("%04d", a[i]); }
    else { printf("%d", a[i]); }
}

std::istream &sjtu::operator>>(std::istream &input, sjtu::int2048 &x) {
  std::string s = "";
  char ch = getchar();
  while (ch == ' ' || ch == '\n' || ch == '\r') ch = getchar();
  while (ch != ' ' && ch != '\n' && ch != '\r')
  {
    s += ch;
    ch = getchar();
  }
  x.read(s);
  return input;
}

std::ostream &sjtu::operator<<(std::ostream &output, const sjtu::int2048 &x)
{
  int2048 tmp = x;
  tmp.print();
  return output;
}

sjtu::int2048 sjtu::int2048::operator+() const
{
  return *this;
}

sjtu::int2048 sjtu::int2048::operator-() const
{
  sjtu::int2048 tmp(*this);
  tmp.sgn *= -1;
  return tmp;
}

sjtu::int2048 &sjtu::int2048::operator=(const sjtu::int2048 &val)
{
  delete [] a;
  len = val.len;
  sgn = val.sgn;
  a = new int [len + 5];
  for (int i = 0; i < len; ++i) a[i] = val.a[i];
  return *this;
}

bool sjtu::operator==(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x.sgn != y.sgn) return false;
  if (x.len != y.len) return false;
  for (int i = 0; i < x.len; ++i)
    if (x.a[i] != y.a[i]) return false;
  return true;
}

bool sjtu::operator!=(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x == y) { return false; }
  else { return true; }
}

bool sjtu::operator<(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x.sgn != y.sgn) return x.sgn < y.sgn;
  if (x.len != y.len) return x.len * x.sgn < y.len * y.sgn;
  for (int i = x.len - 1; i >= 0; --i)
    if (x.a[i] != y.a[i]) return x.a[i] * x.sgn < y.a[i] * y.sgn;
  return false;
}

bool sjtu::operator>(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x.sgn != y.sgn) return x.sgn > y.sgn;
  if (x.len != y.len) return x.len * x.sgn > y.len * y.sgn;
  for (int i = x.len - 1; i >= 0; --i)
    if (x.a[i] != y.a[i]) return x.a[i] * x.sgn > y.a[i] * y.sgn;
  return false;
}

bool sjtu::operator<=(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x > y) { return false; }
  else { return true; }
}

bool sjtu::operator>=(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  if (x < y) { return false; }
  else { return true; }
}

sjtu::int2048 sjtu::UnsignedAdd(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  sjtu::int2048 ans;
  ans.sgn = 1;
  ans.len = std::max(x.len, y.len);
  delete [] ans.a;
  ans.a = new int [ans.len + 5];
  for (int i = 0; i < x.len; ++i) ans.a[i] = x.a[i];
  for (int i = x.len; i < ans.len + 5; ++i) ans.a[i] = 0;
  for (int i = 0; i < y.len; ++i) ans.a[i] += y.a[i];
  for (int i = 1; i < ans.len; ++i)
  {
    ans.a[i] += ans.a[i - 1] / sjtu::int2048::base;
    ans.a[i - 1] %= sjtu::int2048::base;
  }
  if (ans.a[ans.len - 1] >= sjtu::int2048::base)
  {
    ans.a[ans.len] = ans.a[ans.len - 1] / sjtu::int2048::base;
    ans.a[ans.len - 1] %= sjtu::int2048::base;
    ++ans.len;
  }
  return ans;
}

sjtu::int2048 sjtu::UnsignedMinus(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  sjtu::int2048 ans(x);
  ans.sgn = 1;
  int borrow = 0;
  for (int i = 0; i < y.len; ++i)
  {
    ans.a[i] -= borrow, borrow = 0;
    if (ans.a[i] < y.a[i])
    {
      ans.a[i] = ans.a[i] + sjtu::int2048::base - y.a[i];
      borrow = 1;
    }
    else
    {
      ans.a[i] -= y.a[i];
    }
  }
  int cur_digit = y.len;
  while(cur_digit < ans.len && borrow != 0)
  {
    if (ans.a[cur_digit] == 0) { ans.a[cur_digit] = sjtu::int2048::base - 1; }
    else
    {
      --ans.a[cur_digit], borrow = 0;
    }
    ++cur_digit;
  }
  while (ans.a[ans.len - 1] == 0 && ans.len >= 2) --ans.len;
  return ans;
}

sjtu::int2048 &sjtu::int2048::add(const sjtu::int2048 &val)
{
  if (sgn == -1)
  {
    if (val.sgn == -1)
    {
      return *this = -sjtu::UnsignedAdd(*this, val);
    }
    else
    {
      if (-(*this) < val)
      {
        return *this = sjtu::UnsignedMinus(val, *this);
      }
      else
      {
        return *this = -sjtu::UnsignedMinus(*this, val);
      }
    }
  }
  else
  {
    if (val.sgn == 1)
    {
      return *this = sjtu::UnsignedAdd(*this, val);
    }
    else
    {
      if (*this > -val)
      {
        return *this = sjtu::UnsignedMinus(*this, val);
      }
      else
      {
        return *this = -sjtu::UnsignedMinus(val, *this);
      }
    }
  }
}

sjtu::int2048 sjtu::add(sjtu::int2048 x, const sjtu::int2048 &y)
{
  return x.add(y);
}

sjtu::int2048 &sjtu::int2048::minus(const sjtu::int2048 &val)
{
  add(-val);
  return *this;
}

sjtu::int2048 sjtu::minus(sjtu::int2048 x, const sjtu::int2048 &y)
{
  return x.minus(y);
}

sjtu::int2048 &sjtu::int2048::operator+=(const sjtu::int2048 &val)
{
  add(val);
  return *this;
}

sjtu::int2048 sjtu::operator+(sjtu::int2048 x, const sjtu::int2048 &y)
{
  x.add(y);
  return x;
}

sjtu::int2048 &sjtu::int2048::operator-=(const sjtu::int2048 &x)
{
  minus(x);
  return *this;
}

sjtu::int2048 sjtu::operator-(sjtu::int2048 x, const sjtu::int2048 &y)
{
  x.minus(y);
  return x;
}

sjtu::int2048 &sjtu::int2048::operator*=(const sjtu::int2048 &val)
{
  sjtu::polynomial x(*this), y(val);
  int sgn_tmp = sgn * val.sgn;
  *this = x.Multiply(y).ToInteger();
  sgn = sgn_tmp;
  return *this;
}

sjtu::int2048 sjtu::operator*(sjtu::int2048 x, const sjtu::int2048 &y)
{
  x *= y;
  return x;
}

int main()
{
}