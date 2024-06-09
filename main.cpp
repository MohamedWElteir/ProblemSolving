#include <bits/stdc++.h>
#include "code.cpp"
using namespace std;

void fileInput() {
    freopen("input.txt", "r", stdin);
}


void fileOutput() {
    freopen("output.txt", "w", stdout);
}

void fileIO() {
    fileInput();
    fileOutput();
}

void ff() {
    az3arina
#ifndef ONLINE_JUDGE
    fileIO();
#endif
}


template <typename T>
T binarySearchOnAnswer(T low, T high, auto isValid) {
    // You have to know how you are going to move after finding the answer
    T result = -1,mid;
    while (low <= high) {
        mid = low + (high - low) / 2;
        if (isValid(mid)) { // T
            result = mid;
            low = mid + 1;
        } else {  // F
            high = mid - 1;
        }
    }
    return result;
}
template <typename T>
T binarySearchOnAnswer(vector<T> arr, auto isValid) { // overloaded for vectors
    T low = 0, high = arr.size();
    T result = 0;
    while (low <= high) {
       T mid = low + (high - low) / 2;
        if (isValid(mid)) { // T
            result = mid;
            low = mid + 1;
        } else {  // F
            high = mid - 1;
        }
    }
    return result;
}


signed main() {
    az3arina
    ff();
    int n,res=0; in n;
    vectorint arr(n);
    read(arr);
    sort(all(arr));
    res = binarySearchOnAnswer(arr,[&arr](int i){return *(arr.end()-1)-arr[i]*2>0;});
    out res;
    die
}