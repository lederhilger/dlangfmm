module induction.types;

import core.stdc.string : memcpy;
import std.complex : Complex;
import std.math : isFinite, sqrt;
import std.array : overlap;

alias CDouble = Complex!double;
enum maxOrder = 32;
enum maxProfundity = 48;

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

void validateOrder(int order) @safe
{
	enforce(order >= 1 && order <= maxOrder, "order must be in [1,...,32]");
}

package(induction) bool overlaps(const CDouble)[] buffer, const(CDouble)[][] others) @safe @nogc pure nothrow
{
	foreach(other; others)
	{
		if (overlap(buffer, other).length)
		{
			return true;
		}
	}
	return false;
}

package(induction) bool buffersOverlap(const(CDouble)[][] buffers) @safe @nogc pure nothrow
{
	foreach (index; 0 .. buffers.length)
	{
		if (overlaps(buffers[index], buffers[index + 1 .. $]))
		{
			return true;
		}
	}
	return false;
}

ulong keyBits(double value) pure nothrow @nogc
{
	if (value == 0.0) {value = 0.0;}
	ulong bits;
	memcpy(&bits, &value, bits.sizeof);
	return bits;
}

struct GeometryKey
{
	ulong x, y, z;

	this(Vectrix offset) pure nothrow @nogc
	{
		x = keyBits(offset.x);
		y = keyBits(offset.y);
		z = keyBits(offset.z);
	}
}

struct TranslationKey
{
	Translation kind;
	GeometryKey geometry;
}

struct CacheStats
{
	size_t entries;
	size_t geometryEntries;
	size_t sharedEntries;
	size_t bytes;
	size_t keyPayloadBytes;
	size_t hits;
	size_t misses;
	size_t evictions;
	size_t budgetBytes;
}

struct Filament
{
	uint start;
	uint end;
	double circulation = 0.0;
}

struct SourceSet
{
	const(Vectrix)[] vertices;
	const(Filament)[] filaments;
}