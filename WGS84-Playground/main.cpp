#include <iostream>
#include <chrono>

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";

    // 1. Capture the start time
    auto start = std::chrono::steady_clock::now();

    // ... Code you want to measure ...

    // 2. Capture the end time
    auto end = std::chrono::steady_clock::now();

    // 3. Calculate the duration and convert to desired units (e.g., milliseconds)
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time elapsed: " << elapsed.count() << " ms" << std::endl;
    
    return EXIT_SUCCESS;
}
