#pragma once
#ifndef ESSENTIALS_H
#define ESSENTIALS_H
template<typename T>
inline void swap(T& a, T& b) {
	T temp = a;
	a = b;
	b = temp;
}
template<typename T>
inline T abs(const T& x) {
	return (x < 0) ? (-x) : x;
}
#endif