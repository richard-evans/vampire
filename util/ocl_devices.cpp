// Prints a list of all OpenCL platforms and devices
// compile with g++ -std=c++11 -lOpenCL or equivalent

#include <CL/cl.hpp>

#include <vector>
#include <iostream>

int main(void)
{
   std::vector<cl::Platform> platforms;
   cl::Platform::get(&platforms);

   if (platforms.size() == 0)
   {
      std::cout << "No platforms found." << std::endl;
      return 0;
   }

   for (unsigned i=0; i<platforms.size(); ++i)
   {
      std::cout << "Platform " << i << ": " << platforms[i].getInfo<CL_PLATFORM_NAME>() << std::endl;

      std::vector<cl::Device> devices;
      platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);

      if (devices.size() == 0)
      {
         std::cout << "  No devices found." << std::endl;
      }
      else
      {
         for (unsigned j=0; j<devices.size(); ++j)
         {
            std::cout << "  Device " << j << ": " << devices[j].getInfo<CL_DEVICE_NAME>() << std::endl;
         }
      }
   }
}
