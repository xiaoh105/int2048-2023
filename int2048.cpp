#include <iostream>
#include <int2048.h>

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
  long long tmp = val;
  while (tmp)
  {
    ++len;
    tmp /= base;
  }
  if (val == 0) len = 1; // 0 has length 1
  a = new int [len];
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
  a = new int [len];
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
  if (sgn == -1) printf("-");
  for (int i = len - 1; i >= 0; --i) printf("%d", a[i]);
}

sjtu::int2048 sjtu::int2048::operator+() const
{
  return *this;
}

sjtu::int2048 sjtu::int2048::operator-() const
{
  sjtu::int2048 tmp(*this);
  return tmp;
}

sjtu::int2048 &sjtu::int2048::operator=(const sjtu::int2048 &val)
{
  delete [] a;
  len = val.len;
  sgn = val.sgn;
  a = new int [len];
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
  sjtu::int2048 ans(x);
  ans.sgn = 1;
  ans.len = std::max(ans.len, y.len);
  for (int i = 0; i < ans.len; ++i)
  {
    ans.a[i] += y.a[i];
    if (i >= 1)
    {
      ans.a[i] += ans.a[i - 1] / sjtu::int2048::base;
      ans.a[i - 1] %= sjtu::int2048::base;
    }
  }
  while (ans.a[ans.len] >= sjtu::int2048::base)
  {
    ans.a[ans.len + 1] = ans.a[ans.len] / sjtu::int2048::base;
    ans.a[ans.len] %= sjtu::int2048::base;
    ++ans.len;
  }
  return ans;
}

sjtu::int2048 sjtu::UnsignedMinus(const sjtu::int2048 &x, const sjtu::int2048 &y)
{
  sjtu::int2048 ans(x);
  ans.sgn = 1;
  int borrow = 0;
  for (int i = 0; i < ans.len; ++i)
  {
    ans.a[i] -= borrow, borrow = 0;
    if (ans.a[i] < y.a[i])
    {
      ans.a[i] = ans.a[i] + sjtu::int2048::base - y.a[i];
      borrow = 1;
    }
    else
    {
      ans.a[i] -= ans.a[i];
    }
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

int main()
{
}