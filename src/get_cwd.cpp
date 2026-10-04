#include <filesystem>
#include <iostream>

int main() {
    std::string cwd = std::filesystem::current_path().parent_path().string();
    std::cout << cwd << std::endl;
  return 0;
}
