#include <bits/stdc++.h>
using namespace std;
#define die return 0;
#define makeVector(Name, type, size) vector<type> Name(size); read(Name)
#define makeArray(Name, type, size) type Name[size]; read(Name)

#define vectorint vector<int>
#define fix cout << fixed;
#define start(n) n.begin()
#define endstart(n) n.end()
#define fprecision(n) fixed << setprecision(n) <<
#define precision(n) setprecision(n)
#define yes cout << "YES\n"
#define no cout << "NO\n"
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
#define nl '\n'
#define in cin >>
#define out cout <<
#define az3arina  ios_base::sync_with_stdio(false);cin.tie(nullptr);cout.tie(nullptr);
#define readn int n; in n;
#define sp ' '
typedef long long ll;
typedef unsigned long long ull;
const int N = 1e6 + 5;
auto log2_floor(long long i) {return i ? __builtin_clzll(1) - __builtin_clzll(i) : 0;}
auto round_to_target(double number, double target){  return round(number / target) * target; } // n!/(n-k)!
auto digits_count(ll number){ return (int)log10(abs(number))+1;}
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
* Time complexity: O(N)
* Auxiliary space: O(1)
*/
bool isPalindrome(const string &s) {
    unsigned int len = s.length();
    for (unsigned int i = 0, j = len - 1; i != len / 2; i++, j--) {
        if (s[i] != s[j])
            return false;

    }
    return true;
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

void prefix_sum(vectorint &arr, int n) {

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
//    vector<int> arr = {0, 0, 0, 0, 0, 0};
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
    stack<char> st;
    unordered_map<char, char> closeToOpen{{')', '('},
                                          {']', '['},
                                          {'}', '{'},};
    for (auto c: str) {
        if (closeToOpen.find(c) != closeToOpen.end()) {
            if (!st.empty() && st.top() == closeToOpen.operator[](c)) {
                st.pop();
            } else
                return false;
        } else
            st.emplace(c);
    }
    return (st.empty());
}

auto isvalidTag(const vector<string> &str) { // needs some modifications
    stack<string> s;
    unordered_map<string, string> closeToOpen{{"EndHeader", "Header"},
                                              {"EndRow",    "Row"},
                                              {"EndCell",   "Cell"},
                                              {"EndTable",  "Table"},};
    for (const auto &k: str) {
        if (closeToOpen.find(k) != closeToOpen.end()) {
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

auto MaxSlidingWindow(vector<int> &arr, int k) {
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

int dfs(vector<string>& grid, int row, int col, vector<vector<bool>>& visited) {
    if (row < 0 || row >= grid.size() || col < 0 || col >= grid[0].size() || grid[row][col] == '#' || visited[row][col])
        return 0;

    visited[row][col] = true;

    int degree = 1;
    degree += dfs(grid, row + 1, col, visited);
    degree += dfs(grid, row - 1, col, visited);
    degree += dfs(grid, row, col + 1, visited);
    degree += dfs(grid, row, col - 1, visited);

    return degree;
}
auto none(){
    int *p = new int(5);
    auto k = 0x250;
    auto c = p;
    // delete c;
    auto x = [&](){
        auto arr = {69,2,3};
        return data(arr);
    };
    kill_dependency(p);
    auto *z = &c;
    out z;
}
auto Ac()
{

}

auto operator + (string x, double y){
    return x.push_back(y);
}
char grid[10][10]; int vis[10][10];
int solve(int sx,int sy,int ex,int ey,int n){

    if (sx <1 or sy <1 or ex >n or ey > n or grid[sx][sy] == '#' or vis[sx][sy]) return -1;
    if (sx == ex and sy == ey) return 1;

    vis[sx][sy]=1;
    int ch1 = solve(sx+1,sy,ex,ey,n);
    int ch2 = solve(sx-1,sy,ex,ey,n);
    int ch3 = solve(sx,sy+1,ex,ey,n);
    int ch4 = solve(sx,sy-1,ex,ey,n);
    vis[sx][sy] = 0;
    return min({ch1,ch2,ch3,ch4});
}
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
