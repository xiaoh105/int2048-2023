#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

namespace sjtu
{
  class int2048;
  class polynomial
  {
  private:
    // 使用NTT计算乘法, root = 6, mod = 7 * (2 ^ 50) + 1
    constexpr static __int128 root = 6;
    constexpr static __int128 inv = 1313549891316395;
    constexpr static __int128 mod = 7881299347898369;
    int len;
    __int128 *a;
    /// NTT 蝶形变换, O(n)
    void ChangeIndex();
    /// 扩展多项式长度, 自动添加前缀0
    void ExtendLen(int);
    /// 对当前多项式进行快速数论变换
    void NTT(int);
    /// 对多项式进行进位处理
    void CalcCarry();
    /// 带模数的快速幂
    friend __int128 pow_mod(__int128, __int128);
    /// 扩展欧几里得算法，用于计算模数
    friend void Extend_GCD(__int128, __int128, __int128 &, __int128 &);
    /// 求取模意义下的乘法逆元
    friend __int128 inverse(__int128);

  public:
    friend class int2048;
    /// 默认构造函数, 默认构造f(x) = 0
    polynomial();
    /// 复制构造函数
    polynomial(const polynomial &);
    /// 移动构造函数
    polynomial(polynomial &&) noexcept;
    /// 利用大整数生成多项式，即将x转换为f(x) = a0 + a1 * x + ...
    explicit polynomial(const int2048 &);
    /// 析构函数
    ~polynomial();
    /// 复制赋值运算
    polynomial &operator=(const polynomial &);
    /// 移动赋值运算
    polynomial &operator=(polynomial &&) noexcept;
    /// 进行多项式乘法
    polynomial &Multiply(polynomial);
    /// 利用多项式生成大整数
    int2048 ToInteger();
  };
  class int2048
  {
  private:
    const static int base = 10000; // 压位的数字
    const static int base_log10 = 4; // 压位的位数
    int len; // 数字长度，不包含前缀0
    int *a; // 储存数据，0-based, 低位存在0
    /// 无符号加法
    friend int2048 UnsignedAdd(const int2048 &, const int2048 &);
    /// 无符号减法
    friend int2048 UnsignedMinus(const int2048 &, const int2048 &);
    /// 在牛顿迭代法之后进行误差调整
    friend void Adjust(const int2048 &, const int2048&, int2048 &, long long);
    /// 运用牛顿迭代法求逆，即求[2^n / x]
    friend int2048 GetInv(const int2048 &, int);
    /// 无符号除法
    int2048 &UnsignedDivide(const int2048 &);
    /// 将当前整数*(base^x)（左移一个block)
    int2048 &operator<<=(int);
    /// 将当前整数/(base^x)（右移一个block)
    int2048 &operator>>=(int);
    /// 左移操作，定义同上
    friend int2048 operator<<(int2048, int val);
    /// 右移操作，定义同上
    friend int2048 operator>>(int2048, int val);

  public:
    friend class polynomial;
    int sgn;
    /// 默认构造函数
    int2048();
    /// 基于long long的构造函数
    int2048(long long);
    /// 基于字符串的构造函数
    explicit int2048(const std::string &);
    /// 复制构造函数
    int2048(const int2048 &);
    /// 移动构造函数
    int2048(int2048 &&) noexcept;
    /// 析构函数
    ~int2048();

    /// 读入一个大整数
    void read(const std::string &);
    /// 输出储存的大整数
    void print() const;

    /// 返回当前数的绝对值
    friend int2048 abs(const int2048 &);

    /// 加上一个大整数
    int2048 &add(const int2048 &);
    /// 返回两个大整数之和
    friend int2048 add(int2048, const int2048 &);

    /// 减去一个大整数
    int2048 &minus(const int2048 &);
    /// 返回两个大整数之差
    friend int2048 minus(int2048, const int2048 &);

    /// 一元的+运算符，即+x
    int2048 operator+() const;
    /// 一元的-运算符，即-x
    int2048 operator-() const;

    /// 复制赋值运算
    int2048 &operator=(const int2048 &);
    /// 移动赋值运算
    int2048 &operator=(int2048 &&) noexcept;

    int2048 &operator+=(const int2048 &);
    friend int2048 operator+(int2048, const int2048 &);

    int2048 &operator-=(const int2048 &);
    friend int2048 operator-(int2048, const int2048 &);

    int2048 &operator*=(const int2048 &);
    friend int2048 operator*(int2048, const int2048 &);
    friend int2048 operator*(int2048, long long);

    int2048 &operator/=(const int2048 &);
    friend int2048 operator/(int2048, const int2048 &);

    int2048 &operator%=(const int2048 &);
    friend int2048 operator%(int2048, const int2048 &);

    friend std::istream &operator>>(std::istream &, int2048 &);
    friend std::ostream &operator<<(std::ostream &, const int2048 &);

    friend bool operator==(const int2048 &, const int2048 &);
    friend bool operator!=(const int2048 &, const int2048 &);
    friend bool operator<(const int2048 &, const int2048 &);
    friend bool operator>(const int2048 &, const int2048 &);
    friend bool operator<=(const int2048 &, const int2048 &);
    friend bool operator>=(const int2048 &, const int2048 &);
  };
} // namespace sjtu

#endif