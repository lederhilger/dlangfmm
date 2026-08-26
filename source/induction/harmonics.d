module induction.harmonics;

import std.math : sqrt;
import induction.types : CDouble, Vectrix, nmIndex;

private void storeNegative(int n, int m, CDouble positive, CDouble[] output) pure nothrow @nogc
{
	double sign = (m & 1) ? -1.0 : 1.0;
	output[nmIndex(n, -m)] = CDouble(sign * positive.re, -sign * positive.im);
}

void regularSolidHarmonics(int order, Vectrix offset, CDouble[] output) nothrow @nogc
{
	immutable size_t count = cast(size_t)(order * order);
	assert(output.length >= count);
	double rSquared = offset.x * offset.x + offset.y * offset.y + offset.z * offset.z;
	if (rSquared < 1.0e-60)
	{
		output[0 .. count] = CDouble(0.0, 0.0);
		output[0] = CDouble(1.0, 0.0);
		return;
	}

	output[0] = CDouble(1.0, 0.0);
	if (order == 1) {return;}

	output[nmIndex(1, 0)] = -offset.x * output[0];
	foreach (n; 2 .. order)
	{
		double denominator = cast(double)(n * n);
		output[nmIndex(n, 0)] = -((2.0 * n - 1.0) * offset.x * output[nmIndex(n - 1, 0)] + rSquared * output[nmIndex(n - 2, 0)]) / denominator;
	}

	immutable CDouble ω = CDouble(-offset.y, offset.z);
	foreach (m; 1 .. order)
	{
		auto sectorial = ω * output[nmIndex(m - 1, m - 1)] / (2.0 * m);
		output[nmIndex(m, m)] = sectorial;
		storeNegative(m, m, sectorial, output);

		if (m + 1 >= order) {continue;}
		auto adjacent = -offset.x * sectorial;
		output[nmIndex(m + 1, m)] = adjacent;
		storeNegative(m + 1, m, adjacent, output);

		foreach (n; m + 2 .. order)
		{
			double denominator = cast(double)((n - m) * (n + m));
			auto value = -((2.0 * n - 1.0) * offset.x * output[nmIndex(n - 1, m)] + rSquared * output[nmIndex(n - 2, m)]) / denominator;
			output[nmIndex(n, m)] = value;
			storeNegative(n, m, value, output);
		}
	}
}

void singularSolidHarmonics(int order, Vectrix offset, CDouble[] output) nothrow @nogc
{
	immutable size_t count = cast(size_t)(order * order);
	assert(output.length >= count);
	double rSquared = offset.x * offset.x + offset.y * offset.y + offset.z * offset.z;
	if (rSquared < 1.0e-60)
	{
		output[0 .. count] = CDouble(0.0, 0.0);
		output[0] = CDouble(1.0, 0.0);
		return;
	}

	double rInverseSquared = 1.0 / rSquared;
	output[0] = CDouble(1.0 / sqrt(rSquared), 0.0);
	if (order == 1) {return;}

	output[nmIndex(1, 0)] = offset.x * rInverseSquared * output[0];
	foreach (n; 2 .. order)
	{
		double lower = cast(double)((n - 1) * (n - 1));
		output[nmIndex(n, 0)] = rInverseSquared * ((2.0 * n - 1.0) * offset.x * output[nmIndex(n - 1, 0)] - lower * output[nmIndex(n - 2, 0)]);
	}

	immutable CDouble ω = CDouble(-offset.y, offset.z);
	foreach (m; 1 .. order)
	{
		auto sectorial = rInverseSquared * (2.0 * m - 1.0) * ω * output[nmIndex(m - 1, m - 1)];
		output[nmIndex(m, m)] = sectorial;
		storeNegative(m, m, sectorial, output);

		if (m + 1 >= order) {continue;}
		auto adjacent = rInverseSquared * (2.0 * m + 1.0) * offset.x * sectorial;
		output[nmIndex(m + 1, m)] = adjacent;
		storeNegative(m + 1, m, adjacent, output);

		foreach (n; m + 2 .. order)
		{
			double lower = cast(double)((n + m - 1) * (n - m - 1));
			auto value = rInverseSquared * ((2.0 * n - 1.0) * offset.x * output[nmIndex(n - 1, m)] - lower * output[nmIndex(n - 2, m)]);
			output[nmIndex(n, m)] = value;
			storeNegative(n, m, value, output);
		}
	}
}