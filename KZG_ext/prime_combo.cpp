#include <set>
#include <cstdio>
int primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
int N_primes = 10;
std::set<long long> s;
int main()
{
	for(int i = 1; i < (1 << N_primes); ++i)
	{
		long long N = 1;
		for(int j = i, pos = 0; j > 0; j /= 2, pos++)
		{
			if(j % 2 == 1)
				N = N * primes[pos];
		}
		s.insert(N);
	}
	float max_gap = 1;
	for(int i = 100; i < (1 << 20); ++i)
	{
		int gap = *s.lower_bound(i);
		if((float)gap / i > max_gap)
		{
			printf("%d %d\n", gap, i);
		}
		max_gap = std::max((float)gap / i, max_gap);
	}
	printf("%f\n", max_gap);
	return 0;
}
