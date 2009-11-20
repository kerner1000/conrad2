package calhoun.util;


public class PrimeUtil {
	/** The first prime number */
	static int INITIAL_SIZE = 1000;
	
	static long maxPrime = -1;
	public static int[] primes;
	
	/** Returns an array of primes that includes all primes up to at least N.  May contain more. */
	public static int[] primesToAtLeastN(int number) {
		if(maxPrime < number) {
			primes = primesLessThan(number);
			maxPrime = number;
		}
		return primes;
	}
	
	/*
	 * Compute prime numbers, after Knuth, Vol 1, Sec 1.3.2, Alg. "P". 
	 * Note that there may be more efficient algorithms for
	 * finding primes.
	 */
	private static int[] primesLessThan(long stop) {
		int[] prime = new int[INITIAL_SIZE];

		prime[0] = 1; // P1 (ignore prime[0])
		prime[1] = 2; // P1 (ignore prime[0])
		int n = 3; // odd candidates
		int j = 1; // numberFound

		boolean isPrime = true; // for 3
		boolean doMore = true;
		do {
			if (isPrime) {
				if (j == INITIAL_SIZE - 1) {
					// Grow array dynamically if needed
					int[] np = new int[INITIAL_SIZE * 2];
					System.arraycopy(prime, 0, np, 0, INITIAL_SIZE);
					INITIAL_SIZE *= 2;
					prime = np;
				}
				prime[++j] = n; // P2
				isPrime = false;
				if(n > stop)
					doMore = false; 
			}
			n += 2; // P4

			for (int k = 2; k <= j && k < INITIAL_SIZE; k++) { // P5, P6, P8
				long q = n / prime[k];
				long r = n % prime[k];
				if (r == 0) {
					break;
				}
				if (q <= prime[k]) { // P7
					isPrime = true;
					break;
				}
			}
		} while (doMore); // P3
		return prime;
	}
}
