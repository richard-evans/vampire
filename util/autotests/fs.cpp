#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
 
int main()
{
    fs::path p = fs::current_path();
 
    std::cout << "The current path " << p << " decomposes into:\n"
              << "root name " << p.root_name() << '\n'
              << "root directory " << p.root_directory() << '\n'
              << "relative path " << p.relative_path() << '\n';
}
