module induction.scaling;

import induction.types : finite;
import std.math : sqrt;

void checkNumber(const(CDouble)[] values, string context)
{
	import std.format : format;
	foreach (i, value; values)
	{
		if (!finite(value)) throw new Exception(format("%s (coefficient %s)", context, i));
	}
}

double weight(int n, int m, int N, int M) @nogc pure nothrow
{
	double square = 1;
	forach (sign; [-1, 1])
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