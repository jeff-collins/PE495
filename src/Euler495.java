import java.math.BigInteger;
import java.util.*;

public class Euler495 {

    public static List<Integer> primes(int max) {
        System.out.println("Setup: Calculating and caching primes up to " + max);

        long t1 = System.currentTimeMillis();

        int sqrt = (int) Math.sqrt(max);
        BitSet bitSet = new BitSet(max + 1);
        int i, j;
        List<Integer> primes = new ArrayList<>();

        i = 2;
        for (j = i * 2; j <= max; j += i) {
            bitSet.set(j);
        }

        for (i = 3; i <= sqrt; i += 2) {
            if (!bitSet.get(i)) {
                int inc = i * 2;
                for (j = i * 3; j <= max; j += inc) {
                    bitSet.set(j);
                }
            }
        }

        for (i = 2; i <= max; i++)
            if (!bitSet.get(i))
                primes.add(i);

        long t2 = System.currentTimeMillis();
        System.out.println("Setup: Calculated " + primes.size() + " primes in " + (t2 - t1) + "ms");
        return primes;
    }

    static int index = 0;

    public static class IntWrapper implements Comparable<IntWrapper> {
        int val;
        int id;

        public IntWrapper(int val) {
            this.val = val;
            this.id = index++;
        }

        public int compareTo(IntWrapper i) {
            return val > i.val ? 1 : val < i.val ? -1 : id > i.id ? 1 : -1;
        }

        public int hashCode() {
            return val;
        }

        public String toString() {
            return String.valueOf(val);
        }
    }

    public static int getNonUniqueSlotCount(int key) {
        State s = stateCache.get(key);
        return s.duplicateWidthIncludingFloor;
    }

    public static int getCombinationMultiplier(int stateKey, int slots, int count) {
        if (count == 0) {
            return 1;
        }
        int coveredSlots = getNonUniqueSlotCount(stateKey);
        int boxes = slots - coveredSlots;

        if (boxes == 0) {
            return 1;
        }

        if (combinationMultiplierCache.containsKey(new Tuple(count, boxes, 0))) {
            return combinationMultiplierCache.get(new Tuple(count, boxes, 0));
        }

        int combinations = combinations(count + boxes - 1, count);

        combinationMultiplierCache.put(new Tuple(count, boxes, 0), combinations);
        return combinations;
    }

    public static int getKey(State s) {
        if (!reverseStateCache.containsKey(s)) {
            System.out.println("bogus");
        }
        return reverseStateCache.get(s);
    }

    static HashMap<Integer, Integer> statePlusAdditionCache = new HashMap<>();
    public static int getKey(int key, int addSequence) {
        int cacheKey = key << 16 | addSequence;
        if (statePlusAdditionCache.containsKey(cacheKey)) {
            return statePlusAdditionCache.get(cacheKey);
        }
        State s = stateCache.get(key);
        State search = new State(s, addSequence);
        int newStateKey = reverseStateCache.get(search);
        statePlusAdditionCache.put(cacheKey, newStateKey);
        return newStateKey;
    }

    public static void addState(int key, State s) {
        System.out.println("Adding state: " + key + "->" + s);
        stateCache.put(key, s);
        reverseStateCache.put(s, key);
    }

    public static class Tuple {
        int a, b, c;

        public Tuple(int a, int b, int c) {
            this.a = a;
            this.b = b;
            this.c = c;
        }

        public boolean equals(Object o) {
            return o instanceof Tuple && ((Tuple) o).a == a && ((Tuple) o).b == b && ((Tuple) o).c == c;
        }

        public String toString() {
            return a + ":" + b + ":" + c;
        }

        public int hashCode() {
            return a + b + c;
        }
    }

    public static class State {
        TreeSet<IntWrapper> sequences;
        HashMap<Integer, Integer> subStates = new HashMap<>();
        // duplicate width without floor
        int duplicateWidth;
        // duplicate width with floor
        int duplicateWidthIncludingFloor;
        // number of duplications without floor
        int duplicationCount;
        int floor;
        // minimum number of units to form this duplicated state without floor
        // e.g. for 2, 2 it would be 6, a 2 and a 4.  For 2, 2, 3 it would be 15.  2, 4, 9.
        int minSize;
        boolean noDuplication;
        // duplication without floor considered
        int singleDuplication;
        int nonFloorState;
        int stateWithoutFloor;

        void recalc() {
            boolean usedFloor = false;
            HashSet<Integer> set = new HashSet<>();
            ArrayList<IntWrapper> nonFloorSequence = new ArrayList<>();
            for (IntWrapper i : sequences) {
                if (!usedFloor && i.val == this.floor) {
                    usedFloor = true;
                } else {
                    set.add(i.val);
                    duplicateWidth += i.val;
                    duplicationCount++;
                    // creating the steps, each new duplication needs to be multiplied
                    // to be a new height.  The problem is that the minimum needs to
                    // go in reverse order, since it's cheaper to make the larger duplications
                    // lowest.  BUG
                    minSize+= duplicationCount * i.val;
                    nonFloorSequence.add(i);
                }
            }

            if (nonFloorSequence.isEmpty()) {
                noDuplication = true;
            }

            if (nonFloorSequence.size() == 1) {
                singleDuplication = nonFloorSequence.get(0).val;
            }

            for (int skip : set) {
                boolean skipped = false;
                TreeSet<IntWrapper> sequence = new TreeSet<>();
                for(int i = 0; i < nonFloorSequence.size(); i++) {
                    if (!skipped && nonFloorSequence.get(i).val == skip) {
                        skipped = true;
                    } else {
                        sequence.add(nonFloorSequence.get(i));
                    }
                    State s = new State(sequence, 0);
                    int stateIndex = getKey(s);
                    this.subStates.put(skip, stateIndex);
                }
            }
            TreeSet<IntWrapper> sequence = new TreeSet<>();
            sequence.addAll(nonFloorSequence);
            State nonFloor = new State(sequence, 0);
            this.stateWithoutFloor = getKey(nonFloor);
            if (this.floor == 1) {
                // can't convert 1 to a floor, just remove it in the nonfloor state
                this.nonFloorState = this.stateWithoutFloor;
            } else {
                // turn off the floor for the nonFloorState
                this.nonFloorState = getKey(new State(this.sequences, 0));
            }
            this.duplicateWidthIncludingFloor = duplicateWidth + floor;
        }

        public State() {
            this.sequences = new TreeSet<>();
            this.floor = 0;
        }

        public State(int sequence, int floor) {
            this.sequences = new TreeSet<>();
            this.sequences.add(new IntWrapper(sequence));
            this.floor = floor;
        }

        public State(Collection<Integer> collection, int floor) {
            this.sequences = new TreeSet<>();
            for (int i: collection ) {
                this.sequences.add(new IntWrapper(i));
            }
            this.floor = floor;
        }

        public State(TreeSet<IntWrapper> collection, int floor) {
            this.sequences = new TreeSet<>();
            this.sequences.addAll(collection);
            this.floor = floor;
        }

        public State(State s, int addSequence) {
            this.sequences = new TreeSet<>();
            this.sequences.addAll(s.sequences);
            this.sequences.add(new IntWrapper(addSequence));
            this.floor = s.floor;
        }

        public State(State original, State newState) {
            this.sequences = new TreeSet<>();
            this.sequences.addAll(original.sequences);
            this.sequences.addAll(newState.sequences);
            if (original.floor != 0) {
                this.floor = original.floor;
            } else {
                this.floor = newState.floor;
            }
        }

        public boolean collectionsEqual(TreeSet<IntWrapper> a, TreeSet<IntWrapper> b) {
            Iterator<IntWrapper> aIter = a.iterator();
            Iterator<IntWrapper> bIter = b.iterator();
            while(aIter.hasNext()) {
                if (!bIter.hasNext() || aIter.next().val != bIter.next().val) {
                    return false;
                }
            }
            if (bIter.hasNext()) {
                return false;
            }
            return true;
        }

        public boolean equals(Object o) {
            return o instanceof State && collectionsEqual(((State) o).sequences, sequences) && ((State) o).floor == floor;
        }

        public int hashCode() {
            return sequences.hashCode();
        }

        public String toString() {
            return this.sequences.toString() + ":(" + floor + ")";
        }
    }

    // from each state given a count the amount of volume to other states
    public static HashMap<Integer, HashMap<Integer, int[]>> transitionCache = new HashMap<>();

    // from full floor state to other states any slot size
    public static HashMap<Tuple, HashMap<Integer, Integer>> stateCountCache = new HashMap<>();

    // from full floor state to other states any slot size
    public static HashMap<Integer, State> stateCache = new HashMap<>();
    // be able to look up the key
    public static HashMap<State, Integer> reverseStateCache = new HashMap<>();

    public static void incMap(HashMap<Integer, Integer> map, int key, int val) {
        int currCount = 0;
        if (map.containsKey(key)) {
            currCount = map.get(key);
        }
        map.put(key, currCount + val);
    }

    public static void decMap(HashMap<Integer, Integer> map, int key, int val) {
        if (!map.containsKey(key)) {
            throw new RuntimeException("Map doesn't contain key: " + key);
        }
        int currCount = map.get(key);
        if (currCount < val) {
            throw new RuntimeException("Trying to decrement key: " + key + " val " + val + " but only has count " + currCount);
        }
        if (currCount == val) {
            map.remove(key);
        } else {
            map.put(key, currCount - val);
        }
    }

    public static HashMap<Tuple, Integer> combinationMultiplierCache = new HashMap<>();

//    the number of ways to distribute n non-distinct objects into k distinct boxes with
//    some possibly left empty is
//
//            (k + n - 1)
//                (n)
//    http://www.cs.columbia.edu/~cs4205/files/CM5.pdf

    public static HashMap<Integer, Integer> factorialFactors(int n) {
        HashMap<Integer, Integer> result = new HashMap<>();
        for (int i = 2; i <= n; i++) {
            HashMap<Integer, Integer> primes = primeFactors(i);
            add(result, primes);
        }
        return result;
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

        // This condition is to handle the case whien n is a prime number
        // greater than 2
        if (n > 2)
            result.put(n, 1);
        return result;
    }

    public static void add(HashMap<Integer, Integer> map1, HashMap<Integer, Integer> map2) {
        for (int prime :
                map2.keySet()) {
            incMap(map1, prime, map2.get(prime));
        }
    }

    public static void subtract(HashMap<Integer, Integer> map1, HashMap<Integer, Integer> map2) {
        for (int prime :
                map2.keySet()) {
            decMap(map1, prime, map2.get(prime));
        }
    }

    public static HashMap<Integer, Integer> factorize(int n) {
        HashMap<Integer, Integer> factors = new HashMap<>();
         // for each potential factor
        int sqrt = (int)Math.sqrt(n);
        for (int factor = 2; factor <= sqrt; factor++) {

            // if factor is a factor of n, repeatedly divide it out
            while (n % factor == 0) {
                incMap(factors, factor, 1);
                n = n / factor;
            }
        }
        if (n > 1) {
            incMap(factors, n, 1);
        }
        return factors;
    }

    static long m(long n) {
        n %= 1000000000;
        if (n < 0) {
            n += 1000000000;
        }
        return n;
    }

    public static int combinations(int n, int r) {
        long value = 1;
        HashMap<Integer, Integer> result = new HashMap<>();
        for (int i = n - r + 1; i <= n; i++) {
            HashMap<Integer, Integer> primes = primeFactors(i);
            add(result, primes);
        }
        for (int i = 2; i <= r; i++) {
            HashMap<Integer, Integer> primes = primeFactors(i);
            subtract(result, primes);
        }
        for (int prime :
                result.keySet()) {
            int count = result.get(prime);
            for (int i = 0; i < count; i++) {
                value = m(value * prime);
            }
        }
        return (int) value;
    }

//    public static long countCombinations(List<Integer> primes, HashMap<Integer, Integer> primeFactors, int slots) {
//        HashMap<Integer, Integer> stateCounts = new HashMap<>();
//        stateCounts.put(slots, 1);
//        for (int prime : primes) {
//            if (primeFactors.containsKey(prime)) {
//                int count = primeFactors.get(prime);
//                HashMap<Integer, Integer> newStateCounts = new HashMap<>();
//                for (int state :
//                        stateCounts.keySet()) {
//                    int origCount = stateCounts.get(state);
//                    int nonUniqueSlotCount = getNonUniqueSlotCount(state);
//                    for (int i = count; i >= 0; i--) {
//                        if (i < count && nonUniqueSlotCount == slots) {
//                            // no unique slots to throw away count
//                            break;
//                        }
//                        HashMap<Integer, Integer> newStates = dropInState(state, i, slots);
//                        long multiplier = getCombinationMultiplier(state, slots, count - i);
//                        for (int key :
//                                newStates.keySet()) {
//                            int newCount = newStates.get(key);
//                            incMap(newStateCounts, key, (int) m(origCount * m(multiplier * newCount)));
//                        }
//                    }
//                }
//                stateCounts = newStateCounts;
//            }
//        }
//        return stateCounts.get(0);
//    }

//    public static HashMap<Integer, Integer> dropInState(int state, int count, int slots) {
//
//        // need to deal with start case everything flat
//
//        if (transitionCache.containsKey(new Tuple(state, count))) {
//            return transitionCache.get(new Tuple(state, count));
//        }
//
//        HashMap<Integer, Integer> stateCounts = new HashMap<>();
//
//        int currSeq = getLastSeq(state);
//
//        if (count != 0 && currSeq != 0 && slots != 0) {
//            for (int i = count; i >= 0; i--) {
//                int remainder = count - i;
//
//                int childKey = removeLastSeq(state);
//
//                int nextSeq = getLastSeq(childKey);
//                if (nextSeq == 0 && remainder > 0) {
//                    // no slots left to use
//                    break;
//                }
//
//                HashMap<Integer, Integer> childStateCounts = dropInState(childKey, remainder, slots - currSeq);
//                HashMap<Integer, Integer> localSeqCounts = calcStatesOnFlatSurface(currSeq, i);
//
//                for (int localKey :
//                        localSeqCounts.keySet()) {
//                    int localCount = localSeqCounts.get(localKey);
//                    if (childStateCounts.size() > 0) {
//                        for (int newChildKey :
//                                childStateCounts.keySet()) {
//                            int childCount = childStateCounts.get(newChildKey);
//                            int newKey = combineKeys(localKey, newChildKey);
//                            incMap(stateCounts, newKey, localCount * childCount);
//                        }
//                    } else {
//                        incMap(stateCounts, localKey, localCount);
//                    }
//                }
//            }
//        }
//        transitionCache.put(new Tuple(state, count), stateCounts);
//        return stateCounts;
//    }
//
    public static HashMap<Integer, Integer> calcStatesOnFlatSurface(int slots, int count) {
        Tuple token = new Tuple(slots, count, count);
        if (stateCountCache.containsKey(token)) {
            return stateCountCache.get(token);
        }
        HashMap<Integer, Integer> stateCounts = calcStatesOnFlatSurfaceTop(slots, count, 0, 0);
        stateCountCache.put(token, stateCounts);
        return stateCounts;
    }

//    Adding state: 0->[]:(0)
//    Adding state: 1->[1]:(1)
//    Adding state: 2->[1, 2]:(1)
//    Adding state: 3->[1, 3]:(1)
//    Adding state: 4->[2]:(0)
//    Adding state: 5->[2, 2]:(0)
//    Adding state: 6->[2]:(2)
//    Adding state: 7->[2, 2]:(2)
//    Adding state: 8->[3]:(0)
//    Adding state: 9->[3]:(3)
//    Adding state: 10->[4]:(0)
//    Adding state: 11->[4]:(4)

    public static HashMap<Integer, Integer> calcStatesOnFlatSurfaceTop(int slots, int count, int prevCount, int currentIterations) {
        calls++;
        if (calls % 10000 == 0) {
            System.out.println("Calls: " + calls);
        }
        HashMap<Integer, Integer> stateCounts = new HashMap<>();

        // terminate current sequence and loop starting from the next level down
        HashMap<Integer, Integer> childStateCounts = calcStatesOnFlatSurfaceIteration(slots, count, prevCount);
        if (currentIterations <= 1) {
            add(stateCounts, childStateCounts);
        } else {
            for (int key :
                    childStateCounts.keySet()) {
                int newKey = getKey(key, currentIterations);
                incMap(stateCounts, newKey, childStateCounts.get(key));
            }
        }

        if (slots == 0) {
            return stateCounts;
        }

        int max = count;
        if (prevCount > 0 && prevCount < count) {
            max = prevCount;
        }

        if (count >= prevCount && prevCount > 0) {
            // we could have a continuation because there's enough to make the same amount again
            // on the next slot (count >= prevCount), let the recursive call terminate it
            int remainingSlots = slots - 1;

            int remainder = count - max;

            HashMap<Integer, Integer> continuationStateCounts =
                    calcStatesOnFlatSurfaceTop(remainingSlots, remainder, max, currentIterations + 1);
            add(stateCounts, continuationStateCounts);
        }

        return stateCounts;
    }

    public static long cacheHits = 0;
    public static long calls = 0;

    public static HashMap<Integer, Integer> calcStatesOnFlatSurfaceIteration(int slots, int count, int prevCount) {

        int max = count;
        if (prevCount > 0 && count >= prevCount) {
            max = prevCount - 1;
        }

        // cache against any value lower than count or count
        Tuple token = new Tuple(slots, count, Math.min(count, max));
        if (stateCountCache.containsKey(token)) {
            cacheHits++;
            if (cacheHits % 10000 == 0) {
                System.out.println("Cache hits: " + cacheHits);
            }
            return stateCountCache.get(token);
        }

        HashMap<Integer, Integer> stateCounts = new HashMap<>();

        if (count == 0) {
            if (slots > 0) {
                stateCounts.put(getKey(new State(slots, slots)), 1);
            } else {
                stateCounts.put(getKey(new State()), 1);
            }
        } else {

            for (int i = max; i >= 0; i--) {

                if (prevCount == 0) {
                    //System.out.println("Starting top iteration " + i);
                }

                int remainder = count - i;
                int remainingSlots = slots - 1;
                if (((double) remainder / (double) remainingSlots) > i) {
                    break;
                }

                HashMap<Integer, Integer> childStateCounts = calcStatesOnFlatSurfaceTop(remainingSlots, remainder, i, 1);
                add(stateCounts, childStateCounts);
            }
        }
        stateCountCache.put(token, stateCounts);
        return stateCounts;
    }

//    public static long uniqueFactorizationsOfFactorial(int factorial, int slots) {
//        List<Integer> primes = primes(factorial);
//        HashMap<Integer, Integer> factors = factorialFactors(factorial);
//        long combos = countCombinations(primes, factors, slots);
//        combos += countCombinations(primes, factors, slots-1);
//        return m(combos);
//    }
//
//    public static long uniqueFactorizationsOfNumber(int number, int slots) {
//        List<Integer> primes = primes(number);
//        HashMap<Integer, Integer> factors = primeFactors(number);
//        long combos = countCombinations(primes, factors, slots);
//        combos += countCombinations(primes, factors, slots-1);
//        return m(combos);
//    }

    static int keyId = 0;
    public static int countStates(int k) {
        int numStates = countStates(k, 0, 0, 0, new Stack<>());
        for (int stateIndex : stateCache.keySet()) {
            State s = stateCache.get(stateIndex);
            s.recalc();
        }
        return numStates;
    }

    static int floors = 0;

    public static int countStates(int slots, int lastSequence, int consumed, int floor, Stack<Integer> stack) {
        int space = slots - consumed;
        int count = 1;
        State s = new State(stack, floor);
        addState(keyId++, s);
        if (floor != 0) {
            floors++;
        }
        int startRange = lastSequence;
        if (lastSequence == 0) {

            stack.push(1);
            count+=countStates(slots, 2, 1, 1, stack);
            stack.pop();

            startRange = 2;
        }
        for(int i = startRange; i <= space; i++) {
            stack.push(i);
            if (floor != 0) {
                count += countStates(slots, i, consumed + i, floor, stack);
            } else if (i==lastSequence) {
                count+=countStates(slots, i, consumed + i, 0, stack);
            } else {
                count+=countStates(slots, i, consumed + i, 0, stack);
                count+=countStates(slots, i, consumed + i, i, stack);
            }
            stack.pop();
        }
        return count;
    }

    public static int computePosition(int position, int remainder[], BigInteger curr, BigInteger prev, int k, Stack<BigInteger> numbers) {
        if (position == 0) {
            if (curr.equals(BigInteger.ONE) || curr.compareTo(prev) != 1) {
                return 0;
            }
            numbers.push(curr);
            int val = nextCol(remainder, curr, k, numbers);
            numbers.pop();
            return val;
        }
        int sum = 0;
        BigInteger newVal = curr;
        int origValue = remainder[position];

        for (int i = remainder[position]; i >= 0; i--) {
            sum += computePosition(position-1, remainder, newVal, prev, k, numbers);

            remainder[position]--;
            newVal = newVal.multiply(BigInteger.valueOf(position+1));
        }
        remainder[position] = origValue;
        return sum;
    }


    public static HashSet<ArrayList> calcPreviousStates(HashSet<ArrayList>prev, int remove, HashMap<Integer, Integer> states) {
        HashSet<ArrayList> newResults = new HashSet<>();
        int iter = 0;
        for (ArrayList<Integer> array : prev) {
            int numArray[] = new int[array.size()];
            for(int i = 0; i < array.size(); i++) {
                numArray[i] = array.get(i);
            }
            for (int i = 0; i < numArray.length; i++) {
                while (numArray[i] % remove == 0) {
                    numArray[i] /= remove;
                }
            }
            Arrays.sort(numArray);
            ArrayList<Integer> list = new ArrayList();
            for(int i = 0; i < numArray.length; i++) {
                list.add(numArray[i]);
            }
            newResults.add(list);
        }
        for (ArrayList<Integer> numArray : newResults) {
            iter++;
            int state = findState(numArray);
            incMap(states, state, 1);
        }
        return newResults;
    }

    static int findState(ArrayList<Integer> nums) {
        int i = 0;
        int floor = 0;
        while(i < nums.size() && nums.get(i) == 1) {
            floor++;
            i++;
        }
        if (i == nums.size()) {
            return getKey(new State(nums.size(), nums.size()));
        }
        int curr = nums.get(i);
        int currSeq = 1;
        i++;
        TreeSet<IntWrapper> seq = new TreeSet<>();
        if (floor > 0) {
            seq.add(new IntWrapper(floor));
        }
        while(i < nums.size()) {
            if (nums.get(i) == curr) {
                currSeq++;
            } else {
                if (currSeq > 1) {
                    seq.add(new IntWrapper(currSeq));
                    currSeq = 1;
                }
            }
            curr = nums.get(i);
            i++;
        }
        if (currSeq > 1) {
            seq.add(new IntWrapper(currSeq));
        }
        State s = new State(seq, floor);
        return getKey(s);
    }

    static HashSet<ArrayList> slowResult = new HashSet<>();

    static long iter = 0;
    public static int nextCol(int remainder[], BigInteger prev, int k, Stack<BigInteger> numbers) {
        if (k == 1) {
            BigInteger lastNumber = BigInteger.ONE;
            for (int i = 2; i <= remainder.length; i++) {
                if (remainder[i-1] > 0) {
                    BigInteger nextVal = BigInteger.valueOf(i).pow(remainder[i-1]);
                    lastNumber = lastNumber.multiply(nextVal);
                }
            }
            if (lastNumber.compareTo(prev) == 1) {
                iter++;
                if (iter % 1000000 == 0) {
                    for (BigInteger num : numbers) {
                        System.out.print(num + " ");
                    }
                    System.out.println(lastNumber + " ");
                }
//                ArrayList numArray = new ArrayList();
//                for (BigInteger num : numbers) {
//                    numArray.add(num.intValue());
//                }
//                numArray.add(lastNumber.intValue());
//                slowResult.add(numArray);
//                for (BigInteger num : numbers
//                     ) {
//                    System.out.print(num + " ");
//                }
//                System.out.println(lastNumber + " ");
                return 1;
            } else {
                return 0;
            }
        }
        return computePosition(remainder.length-1, remainder, BigInteger.ONE, prev, k-1, numbers);
    }

    public static int computeSlowK(HashMap<Integer, Integer> factors, int k) {
        int max = 1;
        for (int key : factors.keySet()) {
            if (key > max) {
                max = key;
            }
        }
        int factorArray[] = new int[max];
        for (int key : factors.keySet()) {
            factorArray[key-1] = factors.get(key);
        }
        return nextCol(factorArray, BigInteger.ONE, k, new Stack<BigInteger>());
    }

    public static int computeSlow(HashMap<Integer, Integer> factors, int k) {
        int total = computeSlowK(factors, k);
        int totalIncludingOne = computeSlowK(factors, k-1);
        System.out.println("[slow] k=" + k + ":" + total + " k=" + (k-1) + ":" + totalIncludingOne);
        return total + totalIncludingOne;
    }

    public static int compute(HashMap<Integer, Integer> factors, int k) {
        int n = factors.get(2);
        int maxDrop = factors.get(3);

//        primeQMax(n, k);
//
//        countStates(k);

        int total = computeK(factors, k, n, maxDrop);
        int totalIncludingOne = computeK(factors, k-1, n, maxDrop);
        System.out.println("[fast] k=" + k + ":" + total + " k=" + (k-1) + ":" + totalIncludingOne);
        return total + totalIncludingOne;
    }

    private static int computeK(HashMap<Integer, Integer> factors, int k, int n, int maxDrop) {
        int totals[] = new int[stateCache.size()];

        for (int stateIndex : stateCache.keySet()) {
            State state = stateCache.get(stateIndex);
            if (state.duplicateWidthIncludingFloor <= k) {
                totals[stateIndex] = stateCombos(n, k, n, stateIndex);
            }
        }

        int maxPrime = 0;
        for (int key : factors.keySet()) {
            if (key > maxPrime) {
                maxPrime = key;
            }
        }

        computeCompositeTransitions(maxDrop, k, 0, 0, 0, false, new Stack<>(), null);

        for (int i = 3; i <= maxPrime; i++) {
            int newTotals[] = new int[stateCache.size()];
            if (factors.containsKey(i)) {
                int primeCount = factors.get(i);
                if (primeCount > 0) {
                    for (int stateIndex : stateCache.keySet()) {
                        State state = stateCache.get(stateIndex);
                        if (state.duplicateWidthIncludingFloor <= k) {
                            if (transitionCache.containsKey(stateIndex)) {
                                HashMap<Integer, int[]> transitions = transitionCache.get(stateIndex);
                                for (int newStateIndex : transitions.keySet()) {
                                    int[] values = transitions.get(newStateIndex);
                                    if (state.duplicateWidthIncludingFloor < k) {
                                        for (int iter = 0; iter <= primeCount; iter++) {
                                            int newTotal = totals[stateIndex] * values[iter];
                                            int mult = getCombinationMultiplier(stateIndex, k, primeCount - iter);
                                            newTotals[newStateIndex] += newTotal * mult;
                                        }
                                    } else {
                                        newTotals[newStateIndex] += totals[stateIndex] * values[primeCount];
                                    }
                                }
                            }
                        }
                    }
//                    System.out.println("totals=" + convertToMap(totals));
                    totals = newTotals;
                }
            }
        }
//        System.out.println("totals=" + convertToMap(totals));
        return totals[0];
    }

    static HashMap<Integer, Integer> convertToMap(int a[]) {
        HashMap<Integer, Integer> map = new HashMap<>();
        for (int i = 0; i < a.length; i++) {
            if (a[i] > 0) {
                map.put(i, a[i]);
            }
        }
        return map;
    }

    static int[] crossValues(int a[], int b[], int n) {
        int newValues[] = new int[n+1];
        for (int i = 0; i <= n; i++) {
            for(int j = 0; j <= i; j++) {
                newValues[i] += a[j] * b[i-j];
            }
        }
        return newValues;
    }

    public static void computeCompositeTransitions(int n, int k, int lastSequence,
                                                   int consumedSlots, int floor, boolean justAddedFloor, Stack<Integer> stack,
                                                   HashMap<Integer, int[]> currentCompositeValues) {
        int space = k - consumedSlots;
        State s = new State(stack, floor);
        int currentStateKey = getKey(s);

        HashMap<Integer, int[]> newCompositeValues = new HashMap<>();
        if (lastSequence == 0) {
            int values[] = new int[n+1];
            values[0] = 1;
            newCompositeValues.put(0, values);
        } else {
            // get out the flat surface states, including converting floor states to non-floor if necessary
            // for each flat surface state, combine with a current state
            // compute the count array and put in the new current values map
            // increment the transition amount for each new state by iterating combo multiplier,
            //          times new value
            // save the transition amount for state to other states
            HashMap<Integer, int[]> flatSurfaceStateValues = new HashMap<>();
            for (int targetFlatState : stateCache.keySet()) {
                State flatSurfaceState = stateCache.get(targetFlatState);
                if (flatSurfaceState.duplicateWidthIncludingFloor <= k) {
                    int flatSurfaceStateArray[] = new int[n + 1];
                    boolean hasValues = false;
                    for (int i = 0; i <= n; i++) {
                        flatSurfaceStateArray[i] = stateCombos(i, lastSequence, i, targetFlatState);
                        if (flatSurfaceStateArray[i] > 0) {
                            hasValues = true;
                        }
                    }
                    if (!hasValues) {
                        continue;
                    }
                    if (!justAddedFloor && flatSurfaceState.floor != 0) {
                        // need to convert this to another state without a floor
                        // since we just added a nonfloor
                        targetFlatState = flatSurfaceState.nonFloorState;
                    }
                    if (!flatSurfaceStateValues.containsKey(targetFlatState)) {
                        flatSurfaceStateValues.put(targetFlatState, flatSurfaceStateArray);
                    } else {
                        int currentArray[] = flatSurfaceStateValues.get(targetFlatState);
                        for (int i = 0; i < currentArray.length; i++) {
                            currentArray[i] += flatSurfaceStateArray[i];
                        }
                    }
                }
            }

            for (int targetFlatState : flatSurfaceStateValues.keySet()) {
                if (lastSequence < consumedSlots) {
                    State flatSurfaceState = stateCache.get(targetFlatState);
                    // now go to all current values and combine
                    for (int compositeStateIndex : currentCompositeValues.keySet()) {
                        State compositeState = stateCache.get(compositeStateIndex);
                        if (flatSurfaceState.duplicateWidthIncludingFloor + compositeState.duplicateWidthIncludingFloor <= k) {
                            if (compositeState.floor == 0 || flatSurfaceState.floor == 0) {
                                // can't both have a floor since that would mean two floors
                                State combo = new State(compositeState, flatSurfaceState);
                                int comboIndex = getKey(combo);
                                int newValues[] = crossValues(flatSurfaceStateValues.get(targetFlatState),
                                        currentCompositeValues.get(compositeStateIndex), n);
                                if (!newCompositeValues.containsKey(comboIndex)) {
                                    newCompositeValues.put(comboIndex, newValues);
                                } else {
                                    int oldValues[] = newCompositeValues.get(comboIndex);
                                    for(int i = 0; i < newValues.length; i++) {
                                        oldValues[i] += newValues[i];
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // just one duplicated sequence
                    if (!newCompositeValues.containsKey(targetFlatState)) {
                        newCompositeValues.put(targetFlatState, flatSurfaceStateValues.get(targetFlatState));
                    } else {
                        int newValues[] = flatSurfaceStateValues.get(targetFlatState);
                        int oldValues[] = newCompositeValues.get(targetFlatState);
                        for(int i = 0; i < newValues.length; i++) {
                            oldValues[i] += newValues[i];
                        }
                    }
                }
            }
        }
        transitionCache.put(currentStateKey, newCompositeValues);

        int startRange = lastSequence;
        if (lastSequence == 0) {

            stack.push(1);
            computeCompositeTransitions(n, k, 1, 1, 1, true, stack, newCompositeValues);
            stack.pop();

            startRange = 2;
        }

        if (lastSequence == 1) {
            // can't do more than 1 sequence of 1 so bump up
            startRange = 2;
        }

        for(int newSequence = startRange; newSequence <= space; newSequence++) {

            stack.push(newSequence);
            if (floor != 0) {
                computeCompositeTransitions(n, k, newSequence, consumedSlots + newSequence, floor, false, stack, newCompositeValues);
            } else if (newSequence==lastSequence) {
                computeCompositeTransitions(n, k, newSequence, consumedSlots + newSequence, 0, false, stack, newCompositeValues);
            } else {
                computeCompositeTransitions(n, k, newSequence, consumedSlots + newSequence, 0, false, stack, newCompositeValues);
                computeCompositeTransitions(n, k, newSequence, consumedSlots + newSequence, newSequence, true, stack, newCompositeValues);
            }
            stack.pop();
        }
    }

    static long[][] pcache = new long[10000+1][30+1];
    public static long p(int n, int k) {
        if (n == 0 && k == 0) {
            return 1;
        }
        if (n < k || k <= 0) {
            return 0;
        }
        if (k == 1) {
            return 1;
        }
        if (pcache[n][k] != 0) {
            return pcache[n][k];
        }

        long total = 0;
        total += p(n-1, k-1);
        total += p(n-k, k);

        pcache[n][k] = total;
        return total;
    }

    public static long q(int n, int k) {
        return p(n - (k*(k-1))/2, k);
    }

    public static int qSingleDuplicatedSequence(int n, int k, int dupSeqLength) {
        int sum = 0;
        int max = n/dupSeqLength;
        for (int i = 1; i <= max; i++) {
            sum += qWithoutX(n - dupSeqLength*i, k - dupSeqLength, i);
        }
        return sum;
    }

    public static long qOneDup(int n, int k) {
        long sum = 0;
        for (int i = 1; i <= n/2; i++) {
            sum += qWithoutX(n - 2*i, k - 2, i);
        }
        return sum;
    }

    public static long qOneTrip(int n, int k) {
        long sum = 0;
        for (int i = 1; i <= n/3; i++) {
            sum += qWithoutX(n - 3*i, k - 3, i);
        }
        return sum;
    }

    public static long qOneQuad(int n, int k) {
        long sum = 0;
        for (int i = 1; i <= n/4; i++) {
            sum += qWithoutX(n - 4*i, k - 4, i);
        }
        return sum;
    }

    public static long qWithoutX(int n, int k, int x) {

        if (n < 0 || k < 0 || x < 0) {
            return 0;
        }

        if (n == 0 && k > 0) {
            return 0;
        }

        if (k == 0) {
            if (n == 0) {
                return 1;
            } else {
                return 0;
            }
        } else if (k == 1) {
            if (n == x) {
                return 0;
            } else {
                return 1;
            }
        }

        if (x > n) {
            return q(n, k);
        }

        long q = q(n, k);
        if (q == 0) {
            return 0;
        }
        long qwx = qWithoutX(n - x, k - 1, x);
        return q - qwx;
    }

    public static long qWithoutSet(int n, int k, Set s) {
        return 0;
    }

    public static long countState(State s, int n, int k) {
        return 0;
    }

    public static long pent(long vals[], int position) {
        int mult = 1;
        int step = 1;
        long sum = 0;
        while (true) {
            int leftPos = position - step * (3 * step - 1) / 2;
            if (leftPos <= 0) {
                break;
            }
            sum += vals[leftPos-1] * mult;
            int rightPos = position - step * (3 * step + 1) / 2;
            if (rightPos <= 0) {
                break;
            }
            sum += vals[rightPos-1] * mult;
            step++;
            mult*=-1;
        }
        return sum;
    }

    public static int stairs(int k) {
        return k * (k+1) / 2;
    }

    static void primeQMax(int n, int k) {
        for (int width = 2; width < k; width++) {
            for (int height = n/k; height <= n/width; height++) {
                qMax(n-width*height, k-width, height-1);
                if (height % 100 == 0) {
                    System.out.println("width = " + width + " height=" + height + " complete");
                }
            }
        }
    }

    static long clippedCount = 0;
    static int qMaxCache[][][] = new int[10000+1][30+1][5000+1];
    public static int qMax(int n, int k, int max) {
        if (n <= 0) {
            return 0;
        }

        if (k == 1) {
            if (n <= max) {
                return 1;
            } else {
                return 0;
            }
        }

        int val = 0;
        int minRemainder = stairs(k-1);

        // what's the highest next strip we can make with our n?
        // if max is greater, then not clipped
        if (max >= n-minRemainder) {
            val = (int)q(n, k);
        } else {
            if (qMaxCache[n][k][max] != 0) {
                return qMaxCache[n][k][max];
            }
            for (int i = max; i >= 0; i--) {
                int maxRemainder = minRemainder + ((i - k) * (k - 1));
                if (n - i > maxRemainder) {
                    break;
                }
                val += qMax(n - i, k - 1, i - 1);
            }
            qMaxCache[n][k][max] = val;
            clippedCount++;
            if (clippedCount % 1000000 == 0) {
                System.out.println("clipped count = " + clippedCount);
            }
        }

        return val;
    }

//    static HashMap<Tuple, HashMap<Integer, Integer>> widthWithMaxDupRange = new HashMap<>();
//    static int compIterations = 0;
//    // this doesn't work because i can't add in the original state values during the
//    // width computation since the context is lost.
//    static void computeTransitions(int n, int k) {
//
//        // set up the list of relevant states
//        ArrayList<Integer> stateList[] = new ArrayList[k+1];
//        for (int i = 1; i <= k; i++) {
//            stateList[k] = new ArrayList<>();
//            for (int stateIndex : stateCache.keySet()) {
//                int count = stateCombos(n, k, n, stateIndex);
//                if (count > 0) {
//                    stateList[k].add(stateIndex);
//                }
//            }
//        }
//
//        for (int width = 1; width <= k; width++) {
//            HashMap<Integer, Integer> transitions = new HashMap<>();
//            for (int currentCount = 1; currentCount <= n; currentCount++) {
//                for (int prevWidth = 1; prevWidth <= width; prevWidth++) {
//                    HashMap<Integer, Integer> localMaxTransitions = new HashMap<>();
//
//                    int addWidth = width - prevWidth;
//
//                    for (int candidateMax = 1; candidateMax <= addWidth; candidateMax++) {
//                        compIterations++;
//                        HashMap<Integer, Integer> existingTransitions =
//                                widthWithMaxDupRange.get(new Tuple(prevWidth, candidateMax, currentCount));
//                        if (existingTransitions != null) {
//                            for (int stateIndex : stateList[addWidth]) {
//                                int count = stateCombos(n, k, n, stateIndex);
//                                if (count > 0) {
//                                    for (int existingStateIndex : existingTransitions.keySet()) {
//                                        int targetState = getKey(existingStateIndex, addWidth);
//                                        incMap(localMaxTransitions, targetState, count * existingTransitions.get(existingStateIndex));
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    widthWithMaxDupRange.put(new Tuple(width, addWidth, currentCount), localMaxTransitions);
//                }
//            }
//        }
//    }
//
    static HashMap<Long, Integer> initialStateCache = new HashMap<>(1000000);
    public static int stateCombos(int n, int k, int max, int stateIndex) {

        State state = stateCache.get(stateIndex);
        if (n == 0) {
            if (state.duplicateWidthIncludingFloor == k && state.floor == k) {
                return 1;
            } else {
                return 0;
            }
        }
        if (state.floor > 0) {
            k -= state.floor;
            state = stateCache.get(state.stateWithoutFloor);
        }

        long key = (long)n << 48 | (long)k << 32 | (long)max << 16 | (long)stateIndex;
        if (initialStateCache.containsKey(key)) {
            return initialStateCache.get(key);
        }

        if (k <= 0) {
            return 0;
        }

        if (k < state.duplicateWidth) {
            return 0;
        }

        int sum = 0;

        if (state.noDuplication) {
            sum = qMax(n, k, max);
        } else if (state.singleDuplication > 0 && max >= n - k + stairs(k-state.duplicateWidth)) {
            sum = qSingleDuplicatedSequence(n, k, state.singleDuplication);
        } else {
            sum = qWithDupSequences(n, k, max, state);
            initialStateCache.put(key, sum);
        }
        return sum;
    }

    private static int qWithDupSequences(int n, int k, int max, State state) {
        int sum = 0;
        for (int leftSeqWidth : state.subStates.keySet()) {
            int subStateIndex = state.subStates.get(leftSeqWidth);
            State subState = stateCache.get(subStateIndex);

            int lowerBoundWidth = leftSeqWidth;
            int upperBoundWidth = k - subState.duplicateWidth;
            for (int currentWidth = lowerBoundWidth; currentWidth <= upperBoundWidth; currentWidth++) {
                int nonDupArmRange = k - currentWidth - subState.duplicateWidth;
                int minArmSize = stairs(nonDupArmRange) + subState.minSize +
                        nonDupArmRange /* width */ * subState.duplicationCount /* height */;
                int minHatSize = stairs(currentWidth - leftSeqWidth);
                for (int currentHeight = n / currentWidth; currentHeight > 0; currentHeight--) {
                    int remainder = n - (currentWidth * currentHeight);
                    if (remainder >= minArmSize + minHatSize) {
                        // i'm suspicious of an off by one error in the block computation below
                        int maxArmSize = stairs(nonDupArmRange) + subState.minSize + nonDupArmRange * subState.duplicateWidth +
                                (k - currentWidth) * (currentHeight - 1 - nonDupArmRange - subState.duplicationCount);
                        for (int hatCount = remainder - minArmSize, armCount = minArmSize;
                                    hatCount >= minHatSize && armCount <= maxArmSize;
                                    hatCount--, armCount++) {
                            long qHat = qMax(hatCount, currentWidth - leftSeqWidth, max - currentHeight);
                            long qArm = stateCombos(armCount, k - currentWidth, currentHeight - 1, subStateIndex);
                            if (hatCount == 0) qHat = 1;
                            if (armCount == 0) qArm = 1;
                            sum += qHat * qArm;
                        }
                    }
                }
            }
        }
        return sum;
    }

    public static void main(String[] args) {


//        HashMap<Integer, Integer> factors = factorialFactors(10);
//        int slowTotal = computeSlow(factors, 4);
//        System.out.println("total=" + slowTotal);
//        int total2 = computeSlow(factors, 3);
//        System.out.println("total2=" + total2);
//        System.exit(0);

        //HashMap<Integer, Integer> factors = primeFactors(864);
        HashMap<Integer, Integer> factors = factorialFactors(14);
        int k = 7;
        countStates(k);
        int slowTotal = computeSlow(factors, k);
//        int slowTotal = computeSlowK(factors, k);
//        HashMap<Integer, Integer> states = new HashMap<>();
//        HashSet<ArrayList> prevValues = calcPreviousStates(slowResult, 11, states);
//        System.out.println("states=" + states);
//        states = new HashMap<>();
//        prevValues = calcPreviousStates(prevValues, 7, states);
//        System.out.println("states=" + states);
//        states = new HashMap<>();
//        prevValues = calcPreviousStates(prevValues, 5, states);
//        System.out.println("states=" + states);
//        states = new HashMap<>();
//        prevValues = calcPreviousStates(prevValues, 3, states);
//        System.out.println("states=" + states);
//        states = new HashMap<>();
//        prevValues = calcPreviousStates(prevValues, 2, states);
//        System.out.println("states=" + states);
        int total = compute(factors, k);
        System.out.println("total=" + total + " slowTotal=" + slowTotal);
        System.exit(0);


        //List<Integer> primes = primes(144);
//        HashMap<Integer, Integer> factors = primeFactors(144);
        //HashMap<Integer, Integer> factors = factorialFactors(10000);

        int n = factors.get(2);
        int maxDrop = factors.get(3);

        primeQMax(n, k);

        countStates(k);
        for (int stateIndex : stateCache.keySet()) {
            stateCombos(n, k, n, stateIndex);
        }

        computeCompositeTransitions(maxDrop, k, 0, 0, 0, false, new Stack<>(), null);

        for (int i = 3; i <= n; i++) {
            if (factors.containsKey(i)) {
                int primeCount = factors.get(i);
            }
        }
        // baseline ready

        System.exit(0);

        long num = 0;
        long x = qOneDup(5, 3);

        for (int i = 1; i <= 20; i++) {
            System.out.println("\n******" + i + "******");
            long qOneDupVals[] = new long[200];
            long qOneTripVals[] = new long[200];
            long qOneQuadVals[] = new long[200];
            long qVals[] = new long[200];
            long pVals[] = new long[200];
            for (int j = 1; j <= 200; j++) {
                long qOneDup = qOneDup(j, i);
                long qOneTrip = qOneTrip(j, i);
                long qOneQuad = qOneQuad(j, i);
                long q = q(j, i);
                long p = p(j, i);
                qOneDupVals[j-1] = qOneDup;
                qOneTripVals[j-1] = qOneTrip;
                qOneQuadVals[j-1] = qOneQuad;
                qVals[j-1] = q;
                pVals[j-1] = p;
                long qOneDupPent = pent(qOneDupVals, j);
                long qOneTripPent = pent(qOneTripVals, j);
                long qOneQuadPent = pent(qOneQuadVals, j);
                long qPent = pent(qVals, j);
                long pPent = pent(pVals, j);
                System.out.println("(" + j + "):1d=" + qOneDup + "," + (qOneDup-qOneDupPent) +
                        " 1t=" + qOneTrip + "," + (qOneTrip-qOneTripPent) +
                        " 1q=" + qOneQuad + "," + (qOneQuad-qOneQuadPent) +
                        " p=" + p + "," + (p-pPent) +
                        " q=" + q + "," + (q-qPent)
                );
            }
        }

        System.exit(0);

        for (int i = 1; i <= 100; i++) {
            num = qWithoutX(100, 10, i);
            System.out.println("q(100,10) without " + i + " = " + num);
            num = q(100, 10);
            System.out.println("q(100,10) = " + num);
            num = p(100, 10);
            System.out.println("p(100,10) = " + num);
        }
        num = q(100, 10);
        System.out.println("q(100,10) = " + num);
        System.exit(0);

        num = p(10000, 30);
        System.out.println("p = " + num);
        num = q(10000, 30);
        System.out.println("q = " + num);
        num = qWithoutX(9998, 28, 2);
        System.out.println("q without x = " + num);
        num = 0;
        for (int i = 0; i < 5000; i++) {
            num += qWithoutX(10000 - i * 2, 28, i);
        }
        System.out.println("single dup map = " + num);

        System.exit(0);

        // 0
        // 1+1
        // 12+1
        // 13+1
        // 2
        // 2+2
        // 22
        // 22+2
        // 3
        // 3+3
        // 4
        // 4+4

        // -- 4 --
        // 3/3
        // 2/2
        // 2,2/2
        // 2,1/1
        // 4/0

        // -- 5 --
        // 4/4
        // 3/3 - 2
        // 3,1/1
        // 2,2/2 - 2
        // 5/0


        // -- 6 --
        // 600000
        // 510000
        // 420000
        // 411000
        // 330000
        // 321000
        // 311100
        // 222000
        // 221100
        // 211110
        // 111111

        // 5/5
        // 4/4 - 2
        // 3,2/3
        // 4,2/4
        // 3/3
        // 3,2/2
        // 3,3/3
        // 2,2,2/2
        // 1,4/1
        // 6/0

        // -- 10 --
        // A000 - 3/3
        // 9100 - 2/2
        // 8200 - 2/2
        // 8110 - 1,2/1
        // 7300 - 2/2
        // 7210 - 1/1
        // 7111 - 3/0
        // 6400 - 2/2
        // 6310 - 1/1
        // 6220 - 1,2/1
        // 6211 - 2/0
        // 5500 - 2,2/2
        // 5410 - 1/1
        // 5320 - 1/1
        // 5311 - 2/0
        // 5221 - 2/0
        // 4420 - 1,2/1
        // 4411 - 2,2/0
        // 4330 - 1,2/1
        // 4321 - 0/0
        // 4222 - 3/0
        // 3331 - 3/0
        // 3322 - 2,2/0


        // 3/3 - 1 *
        // 2/2 - 4 *
        // 1,2/1 - 4 *
        // 1/1 - 4
        // 3/0 - 3 *
        // 2/0 - 3
        // 2,2/2 - 1 *
        // 2,2/0 - 2
        // 0/0 - 1



        countStates(4, 0, 0, 0, new Stack<>());

        HashMap<Integer, Integer> stateCounts;
        stateCounts = calcStatesOnFlatSurface(4, 4);
        System.out.print("{");
        for (int key : stateCounts.keySet()
             ) {
            System.out.print(stateCache.get(key) + "->" + stateCounts.get(key) + ", ");
        }
        System.out.print("}");

        computeCompositeTransitions(4, 4, 0, 0, 0, false, new Stack<>(), null);

        List<Integer> primes = primes(144);
        HashMap<Integer, Integer> factors2 = primeFactors(144);

        int slots = 4;

//        HashMap<Integer, Integer> states = new HashMap<>();
//        states.putAll(stateCounts);
//        primes.remove(0);
//        for (int prime : primes) {
//            if (factors.containsKey(prime)) {
//                int primeCount = factors.get(prime);
//                for (int i = 0; i <= primeCount; i++) {
//                    HashMap<Integer, Integer> newStates = new HashMap<>();
//                    for (int key :
//                            states.keySet()) {
//                        Tuple token = new Tuple(key, i, 0);
//                        HashMap<Integer, Integer> transitionCounts;
//                        if (i > 0) {
//                            transitionCounts = transitionCache.get(token);
//                        } else {
//                            transitionCounts = new HashMap<>();
//                            transitionCounts.put(key, 1);
//                        }
//                        if (transitionCounts != null) {
//                            int count = states.get(key);
//                            for (int transitionStateKey :
//                                    transitionCounts.keySet()) {
//                                int transitionCount = transitionCounts.get(transitionStateKey);
//                                int comboVal = getCombinationMultiplier(key, slots, primeCount - i);
//                                int newCount = count * transitionCount * comboVal;
//                                incMap(newStates, transitionStateKey, newCount);
//                            }
//                        }
//                    }
//                    states = newStates;
//                }
//            }
//        }

//        System.out.println("State 0 = " + states.get(0));


//        stateCounts = calcStatesOnFlatSurface(30, 10000);
//        System.out.println("statecounts = \n" + stateCounts);

//        long combos = uniqueFactorizationsOfNumber(144, 4);
//        System.out.println("144=: " + combos);
//        combos = uniqueFactorizationsOfFactorial(100, 30);
//        System.out.println("100!=: " + combos);
//        combos = uniqueFactorizationsOfFactorial(10000, 30);
//        System.out.println("10000!=: " + combos);
    }
}
