import java.math.BigInteger;
import java.util.*;

// Use Polya Enumeration Theorem
//
// https://en.wikipedia.org/wiki/P%C3%B3lya_enumeration_theorem
// https://en.wikipedia.org/wiki/Cycle_index
//
// https://math.stackexchange.com/questions/1250364/multiplication-partitioning-into-k-distinct-elements
// https://math.stackexchange.com/questions/1251256/in-how-many-ways-can-1000000-be-expressed-as-a-product-of-five-distinct-positi
//



public class PE495 {

    static final int MODULO = 1000000007;
    static final int MAX_SLOTS = 100;

    public static void main(String[] args) {
        int slots = 30;
        int n = 10000;
//        answer: 789107601

//        int slots = 10;
//        int n = 100;
//        answer: 287549200

        Map<Integer, Integer> factors = factorialFactors(n);

        long time1 = System.currentTimeMillis();
        BigInteger result = polyaEnumerationTheoremCount(slots, factors);
        long time2 = System.currentTimeMillis();

//        System.out.println("result for factors[" + factors + "]=" + result);
        System.out.println(result);
        System.out.println("mod 1000000007=" + result.mod(BigInteger.valueOf(MODULO)));
        System.out.println("in " + (time2-time1) + "ms");
    }

    static BigInteger polyaEnumerationTheoremCount(int slots, Map<Integer, Integer> factors) {
        if (slots > MAX_SLOTS) {
            throw new RuntimeException("slots too high");
        }
        boolean hasASingle = false;
        for (int i : factors.values()) {
            if (i == 1) {
                hasASingle = true;
            }
        }
        CycleIndex ci = cycleIndex(slots);
        System.out.println(ci);
        BigInteger result = BigInteger.ZERO;
        BigInteger nFactorial = BigInteger.valueOf(1);
        for (int i = 2; i <= slots; i++) {
            nFactorial = nFactorial.multiply(BigInteger.valueOf(i));
        }
        int counter = 1;
        for (TreeMap<Integer, Integer> term : ci.terms) {
            if (hasASingle && !term.containsKey(1)) {
                // terms without a[1] will always miss x^1, which means 0 results
                // for any prime that has exactly 1 count.  Not all possibilities
                // are like this, so this is a faulty assumption for general cases.
                // For large factorials with at least one prime ^ 1 it works.
                counter++;
                continue;
            }
            BigInteger divisor = getTermDivisor(term, ci);
            BigInteger multiplier = nFactorial.divide(divisor);

            BigInteger coefficientProduct = BigInteger.ONE;
            BigInteger[] coefficients = extractCoefficients(term, factors.get(2));
            for (int prime : factors.keySet()) {
                coefficientProduct = coefficientProduct.multiply(coefficients[factors.get(prime)]);
            }

            BigInteger termResult = multiplier.multiply(coefficientProduct);
            result = result.add(termResult);
            System.out.println(counter + " term[" + printTerm(term) + "]=" + termResult);
            counter++;
        }

        result = result.divide(nFactorial);
        return result;
    }

    // Extract coefficient from a taylor series identified by the encoded "term" parameter
    // which is a hashmap of a[n] terms keyed to exponential values.  This represents
    // (1/(1-x^n)^v)(...)(...) for all the key/value pairs in the hashmap.  This generates
    // a taylor series that has a coefficient of [x^n] A(X) =

    // Use simple compounding to generate the series with multiple terms.
    static BigInteger[] extractCoefficients(TreeMap<Integer, Integer> term, int n) {
        BigInteger[] coefficients = new BigInteger[n+1];
        coefficients[0] = BigInteger.ONE;

        for (int monomial : term.descendingKeySet()) {
            int pow = term.get(monomial);
            for (int i = 0; i < pow; i++) {
                for(int j = monomial; j <= n; j++) {
                    if (coefficients[j] == null) {
                        coefficients[j] = BigInteger.ZERO;
                    }
                    if (coefficients[j-monomial] != null) {
                        coefficients[j] = coefficients[j].add(coefficients[j - monomial]);
                    }
                }
            }
        }

        return coefficients;
    }

    // Polya Enumeration Theorem against the Z(S[n]) symmetic symmetry group for 30 slots, and number of objects
    // enumerated as the number of distinct primes in the inital "n" value.
    //
    // https://math.stackexchange.com/questions/1250364/multiplication-partitioning-into-k-distinct-elements
    // https://math.stackexchange.com/questions/1251256/in-how-many-ways-can-1000000-be-expressed-as-a-product-of-five-distinct-positi
    //
    // the pattern is that the divisor of an encoded term (a1^n1 x a2^n2 x ...) will be
    //
    // 1) for each a[] monomial, multiply (index ^ power) * power!
    // 2) a parity for inclusion/exclusion negative for terms
    private static BigInteger getTermDivisor(Map<Integer, Integer> term, CycleIndex ci) {
        BigInteger divisor = BigInteger.ONE;
        int totalPow = 0;
        for (int i : term.keySet()) {
            int pow = term.get(i);
            BigInteger p = BigInteger.valueOf(i).pow(pow);
            divisor = divisor.multiply(p);
            for (int j = 2; j <= pow; j++) {
                divisor = divisor.multiply(BigInteger.valueOf(j));
            }

            totalPow += pow;
        }
        if (ci.n % 2 != totalPow % 2) {
            divisor = divisor.negate();
        }
        return divisor;
    }

    // placeholder for a Polya Enumeration Theorem Cycle Index
    static class CycleIndex {
        HashSet<TreeMap<Integer, Integer>> terms = new HashSet<>();
        int n;

        public String toString() {
            StringBuilder b = new StringBuilder();
            boolean first = true;
            for (Map<Integer, Integer> term : terms) {
                if (first) {
                    first = false;
                } else {
                    b.append(" + ");
                }
                b.append("[");
                BigInteger nFactorial = BigInteger.valueOf(1);
                for (int i = 2; i <= n; i++) {
                    nFactorial = nFactorial.multiply(BigInteger.valueOf(i));
                }
                b.append(nFactorial.divide(getTermDivisor(term, this)));
                b.append("]");
                b.append(printTerm(term));
            }
            return b.toString();
        }
    }

    private static CycleIndex cycleIndexCache[] = new CycleIndex[MAX_SLOTS];
    static {
        cycleIndexCache[0] = new CycleIndex();
        cycleIndexCache[0].terms.add(new TreeMap<>());
    }
    // Generate the possible members of the symmetry group S[n]
    // https://en.wikipedia.org/wiki/Cycle_index
    // This is handled by the recurrence relationship:
    //      Z(S[n])=(1/n)sum(i = 1..n, (-1)^i * a[i] * Z(S[n−i])) where Z(S[0])=1
    // Also, inclusion/exclusion for uniqueness of factors is included with the -1 ^ i.
    // Note that the coefficient generation is not handled here, and is deferred
    // because there's an easier pattern to follow that doesn't require fractions and
    // only uses integers in #getTermDivisor.
    static CycleIndex cycleIndex(int n) {
        if (n == 0) {
            return cycleIndexCache[0];
        }
        if (cycleIndexCache[n] != null) {
            return cycleIndexCache[n];
        }
        CycleIndex ci = new CycleIndex();
        ci.n = n;
        for (int i = 1; i <= n; i++) {
            CycleIndex subCi = cycleIndex(n-i);
            for (Map<Integer, Integer> term : subCi.terms) {
                TreeMap<Integer, Integer> newterm = new TreeMap<>();
                newterm.putAll(term);
                if (newterm.containsKey(i)) {
                    newterm.put(i, newterm.get(i) + 1);
                } else {
                    newterm.put(i, 1);
                }
                ci.terms.add(newterm);
            }
        }
        cycleIndexCache[n] = ci;
        return ci;
    }

    public static HashMap<Integer, Integer> primeFactors(int n) {
        HashMap<Integer, Integer> result = new HashMap<>();
        int twos = 0;
        while (n % 2 == 0) {
            n = n / 2;
            twos++;
        }
        if (twos > 0) {
            result.put(2, twos);
        }

        // n must be odd at this point.  So we can skip one element (Note i = i +2)
        for (int i = 3; i <= Math.sqrt(n); i = i + 2) {
            int count = 0;
            // While i divides n, print i and divide n
            while (n % i == 0) {
                n = n / i;
                count++;
            }
            if (count > 0) {
                result.put(i, count);
            }
        }

        if (n > 2) {
            result.put(n, 1);
        }
        return result;
    }

    public static TreeMap<Integer, Integer> factorialFactors(int n) {
        TreeMap<Integer, Integer> result = new TreeMap<>();
        for (int i = 2; i <= n; i++) {
            HashMap<Integer, Integer> primes = primeFactors(i);
            for (int num : primes.keySet()) {
                int currCount = 0;
                if (result.containsKey(num)) {
                    currCount = result.get(num);
                }
                result.put(num, currCount + primes.get(num));
            }
        }
        return result;
    }

    static String printTerm(Map<Integer, Integer> term) {
        StringBuilder builder = new StringBuilder();
        ArrayList<Integer> termKeys = new ArrayList();
        termKeys.addAll(term.keySet());
        Collections.sort(termKeys);
        for (int monomial : termKeys) {
            int pow = term.get(monomial);
            builder.append("(");
            builder.append("a").append(monomial);
            if (term.get(monomial) > 1) {
                builder.append("^").append(pow);
            }
            builder.append(")");
        }
        return builder.toString();
    }
}
