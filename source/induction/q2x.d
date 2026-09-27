module induction.q2x;

import std.complex : conj;
import std.math : PI, isFinite;
import core.exception : onRangeError;
import std.exception : enforce;
import induction.scaling : weight, checkNumber, checkScale;
import induction.types : CDouble, Vectrix, nmIndex, sub, norm, finite, buffersOverlap, maxOrder;

@safe:

struct scratchQ2X
{
	CDouble[] qPrevious;
	CDouble[] qCurrent;
	CDouble[] pPrevious;
	CDouble[] pCurrent;

	this(int order)
	{
		qPrevious = new CDouble[order + 1];
		qCurrent = new CDouble[order + 1];
		pPrevious = new CDouble[order + 1];
		pCurrent = new CDouble[order + 1];
	}
}

void filamentQ2X(Vectrix x1, Vectrix x2, Vectrix center, int order, CDouble[] output, ref scratchQ2X scratch, double expansionScale = 1)
{
	if (order < 1 || order > maxOrder || output.length < cast(size_t)(order * order) || scratch.qPrevious.length < cast(size_t)(order + 1) || scratch.qCurrent.length < cast(size_t)(order + 1) || scratch.pPrevious.length < cast(size_t)(order + 1) || scratch.pCurrent.length < cast(size_t)(order + 1))
	{
		onRangeError();
	}
	enforce(finite(x1) && finite(x2) && finite(center), "filaments must be finite");

	checkScale(expansionScale);
	const(CDouble)[][5] buffers = [output, scratch.qPrevious, scratch.qCurrent, scratch.pPrevious, scratch.pCurrent];
	if (buffersOverlap(buffers))
	{
		onRangeError();
	}
	
	auto r0 = sub(x1, center);
	auto ru = sub(x2, x1);
	double coefficient = norm(ru) / (4.0 * PI);
	r0 = Vectrix(r0.x / expansionScale, r0.y / expansionScale, r0.z / expansionScale);
	ru = Vectrix(ru.x / expansionScale, ru.y / expansionScale, ru.z / expansionScale);
	CDouble ξ0 = .5 * CDouble(r0.z, r0.y);
	CDouble ξu = .5 * CDouble(ru.z, ru.y);
	CDouble η0 = .5 * CDouble(r0.z, -r0.y);
	CDouble ηu = .5 * CDouble(ru.z, -ru.y);
	double axial0 = r0.x;
	double axialu = ru.x;
	CDouble ξ = ξ0 + ξu;
	CDouble η = η0 + ηu;
	double axial = axial0 + axialu;
	immutable CDouble imaginary = CDouble(0, 1);

	scratch.qPrevious[0] = CDouble(1, 0);
	scratch.pPrevious[0] = CDouble(1, 0);
	enforce(finite(r0) && finite(ru) && isFinite(coefficient), "geometry exceeds double");
	output[0] = CDouble(coefficient, 0);

	foreach (n; 1 .. order)
	{
		auto qPrevious = scratch.qPrevious;
		auto qCurrent = scratch.qCurrent;
		auto pPrevious = scratch.pPrevious;
		auto pCurrent = scratch.pCurrent;
		if (n == 1)
		{
			qCurrent[0] = -axial * qPrevious[0];
		}
		else
		{
			qCurrent[0] = weight(n, 0, n - 1, 1) * imaginary * η * qPrevious[1] - imaginary * ξ * conj(qPrevious[1]) - n * axial * qPrevious[0];
			qCurrent[0] /= n;
		}
		if (n > 1)
		{
			foreach (m; 1 .. n)
			{
				qCurrent[m] = weight(n, m, n - 1, m - 1) * imaginary * ξ * qPrevious[m - 1] - weight(n, m, n - 1, m) * axial * qPrevious[m];
				if (m < n - 1)
				   qCurrent[m] += weight(n, m, n - 1, m + 1) * imaginary * η * qPrevious[m + 1];
				qCurrent[m] /= n;
			}
		}
		qCurrent[n] = weight(n, n, n - 1, n - 1) * imaginary * ξ * qPrevious[n - 1] / n;

		if (n == 1)
		{
			pCurrent[0] = -axial0 * pPrevious[0] + qCurrent[0];
		}
		else
		{
			pCurrent[0] = weight(n, 0, n - 1, 1) * (imaginary * η0 * pPrevious[1] - imaginary * ξ0 * conj(pPrevious[1])) - n * axial0 * pPrevious[0] + qCurrent[0];
		}
		pCurrent[0] /= n + 1;
		if (n > 1)
		{
			foreach (m; 1 .. n)
			{
				pCurrent[m] = weight(n, m, n - 1, m - 1) * imaginary * ξ0 * pPrevious[m - 1] - weight(n, m, n - 1, m) * axial0 * pPrevious[m] + qCurrent[m];
				if (m < n - 1)
				{
					pCurrent[m] += weight(n, m, n - 1, m + 1) * imaginary * η0 * pPrevious[m + 1];
				}
				pCurrent[m] /= n + 1;
			}
		}
		pCurrent[n] = (weight(n, n, n - 1, n - 1) * imaginary * ξ0 * pPrevious[n - 1] + qCurrent[n]) / (n + 1);

		size_t centerIndex = nmIndex(n, 0);
		double signedCoefficient = (n & 1) ? -coefficient : coefficient;
		output[centerIndex] = signedCoefficient * pCurrent[0];
		foreach (m; 1 .. n + 1)
		{
			output[centerIndex - m] = signedCoefficient * pCurrent[m];
			double alternating = (m & 1) ? -1.0 : 1.0;
			output[centerIndex + m] = signedCoefficient * alternating * conj(pCurrent[m]);
		}

		scratch.qPrevious = qCurrent;
		scratch.qCurrent = qPrevious;
		scratch.pPrevious = pCurrent;
		scratch.pCurrent = pPrevious;
	}
	checkNumber(output[0 .. cast(size_t)(order * order)], "coefficient exceeds double");
}