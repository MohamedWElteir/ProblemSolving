#include <bits/stdc++.h>

#define die return 0;
#define makeVector(Name, type, size) vector<type> Name(size); read(Name)
#define fix cout << fixed;
#define start(n) n.begin()
#define endstart(n) n.end()
#define fprecision(n) fixed << setprecision(n) <<
#define precision(n) setprecision(n)
#define yes cout << "YES\n";
#define no cout << "NO\n";
#define all(x) x.begin(),x.end()
#define allarray(x) x,x+n
#define t int t; cin >>t; while(t--)
#define for1(n) for(int i = 1; i <=n; i++)
#define rof1(n) for(int i = n; i>=1; i--)
#define for0(n) for (int i = 0; i < n; ++i)
#define rof0(n) for(int i = n; i > 0; i--)
#define forauto(x) for(auto i : x)
#define read(arr) for(auto &z : arr){cin >>z;}
#define printarray(arr) for(auto &z : arr){cout <<z << sp;}
#define nl "\n"
#define in cin >>
#define out cout <<
#define Faaaaaaaaaaaaaasstttt  ios_base::sync_with_stdio(false);cin.tie(nullptr);cout.tie(nullptr);
#define sp " "
typedef long long int lli;
typedef long long ll;
typedef unsigned long long ull;
using namespace std;

// Time complexity: O(Sqrt(n))
// Space complexity: O(1)
bool isPrime(ll number) {
    if (number <= 1)
        return false;

    // 2 and 3 are prime numbers
    if (number <= 3)
        return true;

    // Check if the number is divisible by 2 or 3
    if (number % 2 == 0 || number % 3 == 0)
        return false;

    // Check for prime numbers of the form 6k ± 1
    for (int i = 5; i * i <= number; i += 6) {
        if (number % i == 0 || number % (i + 2) == 0)
            return false;
    }

    return true;
}


void fileInput() {
    freopen("in.txt", "r", stdin);
}


void fileOutput() {
    freopen("out.txt", "w", stdout);
}

void fileIO() {
    fileInput();
    fileOutput();
}

void ff() {
    Faaaaaaaaaaaaaasstttt
#ifndef ONLINE_JUDGE
    fileIO();
#endif
}

/* converts the number to binary and removes leading zeroes
* Time Complexity: O(1)
* Auxiliary Space: O(1)
* Reference: GeeksForGeeks
*/
string binaryRepresentation(unsigned n) {
    return bitset<32>(n).to_string().substr(32 - log2(n));
}

/* Function that returns the sum of a range (query) from a "1 BASED" prefix sum array.
* Time complexity: O(1)
* Auxiliary sp: O(1)
*/
ll PrefixSumFromLtoR(ll (&arr)[], ll l, ll r) {
    if (l != 1)
        return (arr[r - 1] - arr[l - 2]);
    else
        return arr[r - 1];
}

// Function that sums from 1 to n in constant time O(1)
ll sumFromOneToN(ll n) {
    return (n * (n + 1)) / 2;
}

/*
 * Function that removes all zeros from real number and return the number without zeros
 * Time complexity: O(N)
 * Auxiliary sp: O(1)
 */
ll removeZeros(ll n) {
    int i = 0;
    string str = to_string(n);
    string res;
    while (i != str.length()) {
        if (str[i] != '0') {
            res += str[i];
        }
        i++;
    }
    return stoll(res);
}

/*
* Time complexity: O(N)
* Auxiliary sp: O(1)
*/
bool isPalindrome(const string &s) {
    unsigned int len = s.length();
    for (unsigned int i = 0, j = len - 1; i != len / 2; i++, j--) {
        if (s[i] != s[j])
            return false;

    }
    return true;
}

// Merges two subarrays of array[].
// First subarray is arr[begin..mid]
// Second subarray is arr[mid+1..end]
void merge(ll array[], ll const left, ll const mid, ll const right) {
    ll const subArrayOne = mid - left + 1;
    ll const subArrayTwo = right - mid;

    // Create temp arrays
    auto *leftArray = new int[subArrayOne], *rightArray = new int[subArrayTwo];

    // Copy data to temp arrays leftArray[] and rightArray[]
    for (auto i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    for (auto j = 0; j < subArrayTwo; j++)
        rightArray[j] = array[mid + 1 + j];

    auto indexOfSubArrayOne = 0, indexOfSubArrayTwo = 0;
    int indexOfMergedArray = left;

    // Merge the temp arrays back into array[left..right]
    while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
        if (leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        } else {
            array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }

    // Copy the remaining elements of
    // left[], if there are any
    while (indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }

    // Copy the remaining elements of
    // right[], if there are any
    while (indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    delete[] leftArray;
    delete[] rightArray;
}

// begin is for left index and end is right index
// of the sub-array of arr to be sorted
void mergeSort(ll array[], ll const begin, ll const end) {
    if (begin >= end)
        return;

    ll mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid);
    mergeSort(array, mid + 1, end);
    merge(array, begin, mid, end);
}

int minMovesToMakeStudentsHappy(int n, vector<int> &chairs) {
    int moves = 0;
    for (int i = 1; i <= n; ++i) {
        if (chairs[i - 1] != i) {
            // Find the student with chair number i
            int j = find(chairs.begin(), chairs.end(), i) - chairs.begin();

            // Swap chairs of students i and j
            swap(chairs[i - 1], chairs[j]);

            ++moves;
        }
    }
    return moves;
}

// Find max sub-array sum
ll maxsubarray(const ll arr[], ll n) {
    ll localmax = 0, globalmax = LLONG_MIN;
    for (int i = 0; i < n; ++i) {
        localmax = max(arr[i], arr[i] + localmax);
        if (localmax > globalmax)
            globalmax = localmax;
    }
    return globalmax;
}

void prefix_sum(vector<int> arr, int n) {

    for1(n) {
        arr[i] += arr[i - 1]; // Calculates the prefix sum
    }

}

vector<int> postfix_sum(const vector<int> &arr, int n) {
    vector<int> postf(n, 0);
    postf[n] = arr[n]; // The last element remains the same
    rof1(n - 1) {
        postf[i] = postf[i + 1] + arr[i]; // Calculates the postfix sum
    }

    return postf;
}

// Partial sum :
//    vector<int> arr = {0,0, 0, 0, 0, 0};
////   add 1 from 1 to 3
//// add 2 from 2 to 4
//// add 3 in all elements
//// output should be: {4 6 6 5 3}
//// 0 0 0 0 0
//// 1 0 0 -1 0
//// 0 2 0 0 -2
//// 3 0 0 0 0 -3
//// 4 2 0 -1 -2 -3
//// 4 6 6 5 3 -3
//    arr[1]++;
//    arr[4]--;
//    arr[2] += 2;
//    arr[5] -= 2;
//    arr[1] += 3;
//    arr[6] -= 3;
//    for (int i = 1; i < 5; ++i) {
//        arr[i]+= arr[i-1];
//    }
//    printarray(arr);

/*       Sliding window
 *    ------------------
 *    compute the first k elements of the array then set it to min/max
 *    then starting from the kth element add the k+1 th and remove the last element  and update till
 *    the end of the array by the equation:
 *      loop( int i = k; i<n ; i++)
 *      current += arr[k] - arr[i-k]
 * */

auto isvalidParantheses(const string &str) {
    stack<char> s;
    unordered_map<char, char> closeToOpen{{')', '('},
                                          {']', '['},
                                          {'}', '{'},};
    for (char c: str) {
        if (closeToOpen.contains(c)) {
            if (!s.empty() && s.top() == closeToOpen.operator[](c)) {
                s.pop();
            } else
                return false;
        } else
            s.emplace(c);
    }
    return (s.empty());
}

auto isvalidTag(const vector<string> &str) { // needs some modifications
    stack<string> s;
    unordered_map<string, string> closeToOpen{{"EndHeader", "Header"},
                                              {"EndRow",    "Row"},
                                              {"EndCell",   "Cell"},
                                              {"EndTable",  "Table"},};
    for (const string &k: str) {
        if (closeToOpen.contains(k)) {
            if (!s.empty() && s.top() == closeToOpen.operator[](k)) {
                s.pop();
            } else
                return false;
        } else
            s.emplace(k);
    }
    return (s.empty());
}

auto MinSlidingWindow(const vector<int> &arr, int k) {
    int n = arr.size();

    // Calculate the sum of the first window of size k
    int current_sum = 0;
    for (int i = 0; i < k; ++i) {
        current_sum += arr[i];
    }

    // Initialize variables to store the minimum sum and its starting index
    int min_sum = current_sum;
    int min_index = 0; // to track the index

    // Slide the window through the array to find the minimum sum
    for (int i = k; i < n; ++i) {
        current_sum += arr[i] - arr[i - k]; // Add the next element and remove the first element
        // in the window
        if (current_sum < min_sum) {
            min_sum = current_sum;
            min_index = i - k + 1; // Update the starting index of the minimum sum window
        }
    }

    return min_sum; // Adding 1 to the index to convert to 1-based indexing
}

auto MaxSlidingWindow(const vector<int> &arr, int k) {
    int n = arr.size();

    // Calculate the sum of the first window of size k
    int current_sum = 0;
    for (int i = 0; i < k; ++i) {
        current_sum += arr[i];
    }

    // Initialize variables to store the maximum sum and its starting index
    int maxSum = current_sum;
    int maxIndex = 0;

    // Slide the window through the array to find the minimum sum
    for (int i = k; i < n; ++i) {
        current_sum += arr[i] - arr[i - k]; // Add the next element and remove the first element in the window
        if (current_sum > maxSum) {
            maxSum = current_sum;
            maxIndex = i - k + 1; // Update the starting index of the maximum sum window
        }
    }

    return maxSum;// Returning the maximum sum of window elements
}

double MaxAverage(const vector<int> &arr, int k) {
    int n = arr.size();

    double currentSum = 0;
    for (int i = 0; i < k; ++i) {
        currentSum += arr[i];
    }
    double maxAvg = (currentSum) / k;


    for (int i = k; i < n; ++i) {
        currentSum += arr[i] - arr[i - k];
        double currentAvg = (currentSum) / k;
        if (currentAvg > maxAvg) {
            maxAvg = currentAvg;
        }
    }

    return maxAvg;
}

// returns the max length of possible balanced parentheses
auto maxBalancedParenthesesLength(const string &s) {
    int maxLength = 0;
    stack<int> stack;
    stack.push(-1);
    for (int i = 0; i < s.length(); i++) {
        if (s[i] == '(') {
            stack.push(i);
        } else {
            stack.pop();
            if (!stack.empty()) {
                maxLength = max(maxLength, i - stack.top());
            } else {
                stack.push(i);
            }
        }
    }
    return maxLength;
}


// find the unpaired element in a vector (1,1,2,2,3) -> 3
auto findUnpaired(const vector<int> &arr) {
    unordered_map<int, int> counts;
    for (int i: arr) {
        counts[i]++;
    }
    for (auto &[k, v]: counts) {
        if (v % 2 == 1) {
            return k;

        }
    }
    return -1;
}

auto twoPointers(vector<int> v, auto target) {
    auto n = v.size();
    int c = 0;
    sort(v.begin(), v.end(), greater());
    for (int i = 0; i < n; ++i) {
        int j = i + 1, temp = 0;
        while (v[i] + v[j] < target && j < n) {
            temp++;
            j++;
            c = max(c, temp);
        }

    }
    return (c == 0) ? 1 : c;
}

int maxBooks(vector<int> a, int target) {
    int n = a.size();
    int left = 0, right = 0;
    int sum = 0;
    int maxBooks = 0;
    while (right < n) {

        sum += a[right];

        while (sum > target) {
            sum -= a[left];
            left++;
        }
        maxBooks = max(maxBooks, right - left + 1);
        right++;
    }
    return maxBooks;
}

string increaseNumberRecursive(string s, int index) {

    for (int i = index; i < s.length(); ++i) {
        s[i] = '0';
    }

    if (index > 0)
        if (s[index - 1] == '9')
            s = increaseNumberRecursive(s, index - 1);
        else
            s[index - 1] = ((int) s[index - 1]) + 1;
    else
        s = "1" + s;

    return s;
}


void B() {
    int lastZeroPosition;
    t {
        string s;
        cin >> s;
        size_t index = s.find_last_of("56789");
        lastZeroPosition = s.length();
        while (index != s.npos) {
            s = increaseNumberRecursive(s, index);
            index = s.find_last_of("56789");
        }
        cout << s << '\n';
    }
}

void C() {
    int tcc;
    cin >> tcc;
    while (tcc--) {
        int n;
        cin >> n;
        int bSize = n * (n - 1) / 2;
        vector<int> b(bSize);
        vector<int> a;
        for (int i = 0; i < bSize; ++i) {
            cin >> b[i];
        }
        sort(b.begin(), b.end());

        int counterOfComparisons = n - 1;
        int index = 0;
        while (index < bSize) {
            a.emplace_back(b[index]);
            index += counterOfComparisons;
            --counterOfComparisons;
        }
        while (a.size() < n) {
            a.push_back(*a.rbegin());
        }
        for (auto num: a) cout << num << ' ';
        cout << '\n';
    }
}

void D() {
/*
 Concept:  au−av≥bu−bv

 which means: au - bu >= av - bv;
 so lets do this preprocessing, each index by the array can be represented by 1 number
 */

    t {
        // input
        int n;
        cin >> n;
        vector<long long> a(n);
        vector<long long> b(n);

        vector<pair<long long, int>>
                nodeWithIndex(n);
        for (int i = 0; i < n; ++i) {
            cin >> a[i];
        }

        long long mxNodePower = LONG_LONG_MIN;

        for (int i = 0; i < n; ++i) {
            cin >> b[i];

            // input preprocessing, represent each node by a single number
            a[i] = a[i] - b[i];
            mxNodePower = max(mxNodePower, a[i]);
            nodeWithIndex[i] = {a[i], i + 1};
        }

        vector<int> result;

        for (auto nodePowerWithIndex: nodeWithIndex) {
            if (nodePowerWithIndex.first == mxNodePower) result.emplace_back(nodePowerWithIndex.second);
        }
        cout << result.size() << '\n';
        for (auto res: result) cout << res << ' ';
        out nl;

    }
}

vector<long long> prefixSumOnVector(vector<long long> &nums, bool isPrefixOneBased = true) {
    vector<long long> prefix(nums.size());
    long long sum = 0;
    prefix[0] = 0;

    for (int i = 0 + isPrefixOneBased; i < nums.size(); ++i) {
        sum += nums[i];
        prefix[i] = sum;
    }
    return prefix;
}

vector<ll> partialSumOneBased(vector<pair<ll, ll>>

                              ranges,
                              int size,
                              bool partialOnly = false
) {
//    int size = ranges.size();
    vector<ll> partial(size + 2, 0), prefix(size + 2, 0);

    ll mn = ranges[0].first, mx = ranges[0].second;

    for (
        pair<ll, ll> range
            : ranges) {
        partial[range.first]++;
        partial[range.second + 1]--;
        mn = min(range.first, mn);
        mx = max(range.second, mx);
    }

    if (partialOnly)
        return
                partial;

    return
            prefixSumOnVector(partial);

}

void E() {
    t {
        int n, k;
        cin >> n;
        vector<int> vec(n);
        vector<pair<long long, long long>>
                prefixArrays = {};

        for (int i = 0; i < n; ++i) {
            cin >> vec[i];
        }


        for (int i = 0; i < n; ++i) {
            int mn, mx, biggestMx;
            for (int j = 0; j < n; ++j) {
                mn = min(vec[i], vec[j]);
                mx = max(vec[i], vec[j]);
                biggestMx = max(biggestMx, mx);

                prefixArrays.push_back({mn, mx});
            }
            vector<long long> result = partialSumOneBased(prefixArrays, biggestMx);
            prefixArrays.clear();

            long long sum = 0;
            for (auto x: result)
                sum += x;
            cout << sum << ' ';
        }
        cout << '\n';
    }
}

auto maximizeProfit(int n, int w, vector<int> &weights) {
    vector<vector<int>> rel(n + 1, vector<int>(w + 1, 0));

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= w; ++j) {
            if (weights[i - 1] <= j) {
                rel[i][j] = max(rel[i - 1][j], rel[i - 1][j - weights[i - 1]] + weights[i - 1]);
            } else {
                rel[i][j] = rel[i - 1][j];
            }
        }
    }

    return rel[n][w];
}

auto Cleopatra(int n) {
    vector<int> primes;
    for (int num = 2; num <= n; ++num) {
        if (isPrime(num)) {
            primes.push_back(num);
        }
    }
    int count = 0;
    for (int prime: primes) {
        if (isPrime(prime + 1)) {
            count++;
        }
    }

    return count;
}

void solve() {
    int n, c = 0;
    in n;
    string s;
    in s;
    while (next_permutation(s.begin() + 1, s.end())) {
        c++;
    }
    out c;
}


void removezeros(vector<int> &nums) {
    int nonzero = 0;
    for (int i = 0; i < nums.size(); ++i) {
        if (nums[i] != 0) {
            swap(nums[nonzero++], nums[i]);
        }
    }
}

void fence() {

    int n;
    in n;

    vector<int> a(n + 1);
    int mx = 0;
    for (int i = 1; i <= n; i++) {
        cin >> a[i];
        mx = max(mx, a[i]);
    }

    if (mx > n) {
        no
        return;
    }

    vector<int> b(n + 1);

    for (int i = 1; i <= n; i++) {
        b[a[i]]++;
    }

    for (int i = n - 1; i > 0; i--) {
        b[i] += b[i + 1];
    }

    for (int i = 1; i <= n; i++) {
        if (a[i] != b[i]) {
            no
            return;
        }
    }

    yes

}

void A() {
    t {
        int x, y, n;
        cin >> x >> y >> n;
        int d = y - x;
        if (d <= n - 1) {
            out -1 << nl;
            continue;
        }
        vector<int> a(n);
        a[0] = x;
        a[n - 1] = y;
        for (int i = 1; i < n - 1; i++) {
            a[i] = a[i - 1] + d / (n - i + 1);
            d -= d / (n - i - 1);
        }
        printarray(a);
        out nl;
    }
}

int findGCDexcludingthekthelement(vector<int> arr, int n, int k) {
    int result = 0;
    for (int i = 0; i < n; i++) {
        if (i + 1 == k) { continue; }
        result = gcd(arr[i], result);


    }
    return result;
}

double_t fact(int n) {
    if (n == 0 || n == 1)
        return 1;
    return n * fact(n - 1);
}

double_t findCatalan(int n) {
    return fact(2 * n) / (fact(n + 1) * fact(n));
}

void coinChangeDFS(vector<int> &coins, int amount, int index, vector<int> &current, vector<vector<int>> &result) {
    if (amount == 0) {
        result.push_back(current); // Found a combination that sums up to the amount
        return;
    }

    for (int i = index; i < coins.size(); ++i) {
        if (coins[i] <= amount) {
            current.push_back(coins[i]); // Include the coin
            coinChangeDFS(coins, amount - coins[i], i, current, result); // Explore further with remaining amount
            current.pop_back(); // Backtrack - remove the last added coin for other combinations
        }
    }
}

vector<vector<int>> coinChange(vector<int> &coins, int amount) {
    vector<vector<int>> result;
    vector<int> current;
    coinChangeDFS(coins, amount, 0, current, result);
    return result;
}
auto fab(ll n){
    map< int , int > memo;
    if (n==1) return 1;
    else if (memo[n]) return memo[n];
    else return fab(n-1) + fab(n-2);
}

int main() {
    Faaaaaaaaaaaaaasstttt
    ff();
    int n; in n;
    out fab(n);
   die
}

