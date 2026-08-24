module induction.types;

import core.stdc.string : memcpy;
import std.complex : Complex;
import std.math : isFinite, sqrt;

alias CDouble = Complex!double;

enum Formulation : ubyte
{
	standard,
	lambHelmholtz
}

enum Translator : ubyte
{
	standard,
	rcr
}

enum Truncation : ubyte
{
	centroid,
	contained,
	intersects
}

enum Translation : ubyte
{
	m2m,
	m2l,
	l2l
}

struct Vectrix
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	this(double xValue, double yValue, double zValue) pure nothrow @nogc
	{
		x = xValue;
		y = yValue;
		z = zValue;
	}
}

Vectrix add(Vectrix a, Vectrix b) pure nothrow @nogc
{
	return Vectrix(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vectrix sub(Vectrix a, Vectrix b) pure nothrow @nogc
{
	return Vectrix(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vectrix scale(Vectrix a, double value) pure nothrow @nogc
{
	return Vectrix(a.x * value, a.y * value, a.z * value);
}

Vectrix wedge(Vectrix a, Vectrix b) pure nothrow @nogc
{
	return Vectrix(
	       a.y * b.z - a.z * b.y,
	       a.z * b.x - a.x * b.z,
	       a.x * b.y - a.y * b.x
	);
}

double dot(Vectrix a, Vectrix b) pure nothrow @nogc
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double normSquared(Vectrix value) pure nothrow @nogc
{
	return dot(value, value);
}

double norm(Vectrix value) pure nothrow @nogc
{
	return sqrt(normSquared(value));
}

bool finite(Vectrix value) pure nothrow @nogc
{
	return isFinite(value.x) && isFinite(value.y) && isFinite(value.z);
}

size_t nmIndex(int n, int m) pure nothrow @nogc
{
	return cast(size_t)(n * n + n + m);
}