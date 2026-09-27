module induction.scaling;

import induction.types : CDouble, finite;
import std.math : sqrt, isFinite;
import std.exception : enforce;

package(induction):
@safe:

void checkNumber(const(CDouble)[] values, string context)
{
	import std.format : format;
	foreach (i, value; values)
	{
		if (!finite(value)) throw new Exception(format("%s (coefficient %s)", context, i));
	}
}

void checkScale(double value)
{
	enforce(isFinite(value) && value > 0, "scale must be positive and finite");
}

double weight(int n, int m, int N, int M) @nogc pure nothrow
{
	double square = 1;
	foreach (sign; [-1, 1])
	{
		int first = n + sign * m;
		int second = N + sign * M;
		foreach (k; second + 1 .. first + 1)
		{
			square *= k;
		}
		foreach (k; first + 1 .. second + 1)
		{
			square /= k;
		}
	}
	return sqrt(square);
}

struct BinomialFactors
{
}

struct RadialFactor
{
}