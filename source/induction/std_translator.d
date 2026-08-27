module induction.std_translator;
import std.exception : enforce;
import induction.harmonics : regularSolidHarmonics, singularSolidHarmonics;
import induction.types : CDouble, Vectrix, nmIndex, Translation, TranslationKey, GeometryKey, CacheStats;

final class TranslationData
{
	ulong generation;
	Translation kind;
	Vectrix offset;
	CDouble[] harmonics;

	this(ulong generationValue, Translation kindValue, Vectrix offsetValue, CDouble[] harmonicValues)
	{
		generation = generationValue;
		kind = kindValue;
		offset = offsetValue;
		harmonics = harmonicValues;
	}
}

private struct SparsePattern
{
	uint[] rowOffsets;
	uint[] terms;
}

private struct CacheRecord
{
	TranslationKey key;
	ulong generation;
	size_t bytes;
}

final class StandardTranslator
{
	private int multipoleOrder;
	private size_t coefficientCount;
	private SparsePattern m2mPattern;
	private SparsePattern l2lPattern;
	private TranslationData[TranslationKey] translationCache;
	private CacheRecord[] evictionRecords;
	private size_t evictionCursor;
	private ulong nextGeneration;
	private CacheStats statistics;

	this(int order, size_t cacheBudgetBytes)
	{
		enforce(order >= 1, "order must be positive");
		enforce(order <= 128, "order too large (max 128)");
		multipoleOrder = order;
		coefficientCount = cast(size_t)(order * order);
		statistics.budgetBytes = cacheBudgetBytes;

		m2mPattern = buildSparsePattern(Translation.m2m);
		l2lPattern = buildSparsePattern(Translation.l2l);
	}

	@property int order() const pure nothrow @nogc { return multipoleOrder; }
	@property size_t persistentBytes() const pure nothrow @nogc
	{
		return (m2mPattern.rowOffsets.length + l2lPattern.rowOffsets.length) * uint.sizeof + (m2mPattern.terms.length + l2lPattern.terms.length) * uint.sizeof;
	}

	private SparsePattern buildSparsePattern(Translation kind)
	{
		auto rowOffsets = new uint[coefficientCount + 1];
		uint[] terms;
		foreach (targetN; 0 .. multipoleOrder)
		{
			foreach (targetM; -targetN .. targetN + 1)
			{
				size_t targetIndex = nmIndex(targetN, targetM);
				rowOffsets[targetIndex] = cast(uint)terms.length;
				int sourceBegin = kind == Translation.m2m ? 0 : targetN;
				int sourceEnd = kind == Translation.m2m ? targetN + 1 : multipoleOrder;
				foreach (sourceN; sourceBegin .. sourceEnd)
				{
					int n = kind == Translation.m2m ? targetN - sourceN : sourceN - targetN;
					int lowerM = -sourceN;
					int upperM = sourceN;
					if (targetM - n > lowerM) {lowerM = targetM - n;}
					if (targetM + n < upperM) {upperM = targetM + n;}
					foreach (sourceM; lowerM .. upperM + 1)
					{
						uint sourceIndex = cast(uint)nmIndex(sourceN, sourceM);
						uint harmonicIndex = cast(uint)nmIndex(n, sourceM - targetM);
						terms ~= sourceIndex | (harmonicIndex << 16);
					}
				}
			}
		}
		rowOffsets[coefficientCount] = cast(uint)terms.length;
		return SparsePattern(rowOffsets, terms);
	}

	private void enforceBudget(TranslationData newest)
	{
		while (statistics.bytes > statistics.budgetBytes && evictionCursor < evictionRecords.length)
		{
			auto record = evictionRecords[evictionCursor++];
			auto found = record.key in translationCache;
			if (found is null || (*found).generation != record.generation) {continue;}
			if (*found is newest)
			{
				--evictionCursor;
				break;
			}
			translationCache.remove(record.key);
			statistics.bytes -= record.bytes;
			--statistics.entries;
			--statistics.geometryEntries;
			statistics.keyPayloadBytes -= GeometryKey.sizeof;
			++statistics.evictions;
		}
		if (evictionCursor > 1024 && evictionCursor * 2 > evictionRecords.length)
		{
			evictionRecords = evictionRecords[evictionCursor .. $].dup;
			evictionCursor = 0;
		}
	}

	TranslationData prepare(Translation kind, Vectrix offset)
	{
		auto key = TranslationKey(kind, GeometryKey(offset));
		auto found = key in translationCache;
		if (found !is null)
		{
			++statistics.hits;
			return *found;
		}
		++statistics.misses;
		auto harmonics = buildHarmonics(kind, offset);
		ulong generation = ++nextGeneration;
		auto cached = new TranslationData(generation, kind, offset, harmonics);
		translationCache[key] = cached;
		immutable size_t bytes = harmonics.length * CDouble.sizeof + GeometryKey.sizeof + __traits(classInstanceSize, TranslationData);
		evictionRecords ~= CacheRecord(key, generation, bytes);
		++statistics.entries;
		++statistics.geometryEntries;
		statistics.keyPayloadBytes += GeometryKey.sizeof;
		statistics.bytes += bytes;
		enforceBudget(cached);
		if (statistics.bytes > statistics.budgetBytes)
		{
			translationCache.remove(key);
			statistics.bytes -= bytes;
			--statistics.entries;
			--statistics.geometryEntries;
			statistics.keyPayloadBytes -= GeometryKey.sizeof;
			++statistics.evictions;
		}
		return cached;
	}

	private CDouble[] buildHarmonics(Translation kind, Vectrix offset)
	{
		int tableOrder = kind == Translation.m2l ? 2 * multipoleOrder - 1 : multipoleOrder;
		auto table = new CDouble[tableOrder * tableOrder];
		if (kind == Translation.m2l) {singularSolidHarmonics(tableOrder, offset, table);}
		else {regularSolidHarmonics(tableOrder, offset, table);}
		return table;
	}

	void apply(int potentialCount)(const TranslationData data, const(CDouble)[] source, CDouble[] target) const nothrow @nogc
	{
		static assert(potentialCount == 2 || potentialCount == 3, "support for std & lh implemented");
		assert(data !is null);
		assert(source.length == cast(size_t)potentialCount * coefficientCount);
		assert(target.length == source.length);
		final switch (data.kind)
		{
			case Translation.m2m:
			{
				applySparse!potentialCount(data, m2mPattern, source, target);
				break;
			}
			case Translation.m2l:
			{
				applyM2L!potentialCount(data, source, target);
				break;
			}
			case Translation.l2l:
			{
				applySparse!potentialCount(data, l2lPattern, source, target);
				break;
			}
		}
	}

	private void applySparse(int potentialCount)(const TranslationData data, const SparsePattern pattern, const(CDouble)[] source, CDouble[] target) const nothrow @nogc
	{
		foreach (targetIndex; 0 .. coefficientCount)
		{
			CDouble[potentialCount] sums = CDouble(0, 0);
			size_t begin = pattern.rowOffsets[targetIndex];
			size_t end = pattern.rowOffsets[targetIndex + 1];
			foreach (termIndex; begin .. end)
			{
				uint packed = pattern.terms[termIndex];
				size_t sourceIndex = packed & 0xffffu;
				auto a = data.harmonics[packed >> 16];
				static foreach (potentialIndex; 0 .. potentialCount)
				{
					sums[potentialIndex] += source[potentialIndex * coefficientCount + sourceIndex] * a;
				}
			}
			static foreach (potentialIndex; 0 .. potentialCount)
			{
				target[potentialIndex * coefficientCount + targetIndex] = sums[potentialIndex];
			}
		}
	}

	private void applyM2L(int potentialCount)(const TranslationData data, const(CDouble[]) source, CDouble[] target) const nothrow @nogc
	{
		size_t targetIndex;
		foreach (targetN; 0 .. multipoleOrder)
		{
			foreach (targetM; -targetN .. targetN + 1)
			{
				CDouble[potentialCount] sums = CDouble(0, 0);
				size_t sourceIndex;
				foreach (sourceN; 0 .. multipoleOrder)
				{
					int harmonicN = targetN + sourceN;
					size_t harmonicIndex = nmIndex(harmonicN, -sourceN - targetM);
					foreach (sourceM; -sourceN .. sourceN + 1)
					{
						auto a = data.harmonics[harmonicIndex++];
						static foreach (potentialIndex; 0 .. potentialCount)
						{
							sums[potentialIndex] += source[potentialIndex * coefficientCount + sourceIndex] * a;
						}
						++sourceIndex;
					}
				}
				static foreach (potentialIndex; 0 .. potentialCount)
				{
					target[potentialIndex * coefficientCount + targetIndex] = sums[potentialIndex];
				}
				++targetIndex;
				
			}
		}
	}

	@property CacheStats stats() const {return statistics;}
}