#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <CL/opencl.h>

#include "errors.h"

#define BUFFER_SIZE 1024
#define SQRSIZE 4096
#define NO_OFFSET 0
#define DIMENSIONS 2

void process_error(cl_int err, const char *function_name) {
  if (function_name != NULL)
    printf("Error calling %s: %s\n", function_name, get_error_string(err));
  else
    printf("Error: %s\n", get_error_string(err));

  exit(EXIT_FAILURE);
}

int main (int argc, char *argv[])
{
  // matmult related variables
  double *h_A, *h_B, *h_C;
  size_t nA, nB, mA, mB;
  size_t sizeA, sizeB, sizeC;
  size_t sizeA_bytes, sizeB_bytes, sizeC_bytes;

  cl_mem d_A, d_B, d_C, d_nA;
  cl_mem_flags flags;
  //

  // opencl launch related variables
  int i, j;
  cl_int err;
  cl_uint num_platforms;
  cl_platform_id *platforms;

  size_t param_value_size_ret;

  cl_uint num_devices;
  cl_device_id *devices;
  cl_device_id device;

  cl_context context;
  cl_command_queue command_queue;

  char *program_filename, *program_buffer;
  FILE *fd_program;
  size_t program_size;
  cl_program program;

  cl_kernel kernel;

  size_t *global_work_size;

  cl_uint work_dim;
  cl_event event;
  //

  // matmult related stuff
  nA = nB = mA = mB = SQRSIZE;
  sizeA = nA * mA;
  sizeB = nB * mB;
  sizeC = nA * mB;
  sizeA_bytes = sizeA * sizeof(double);
  sizeB_bytes = sizeB * sizeof(double);
  sizeC_bytes = sizeC * sizeof(double);

  h_A = (double *) malloc(sizeA_bytes);
  h_B = (double *) malloc(sizeB_bytes);
  h_C = (double *) malloc(sizeC_bytes);

  for (i = 0; i < sizeA; i++) {
    h_A[i] = h_B[i]= 2.;
  }
  // end matmult related stuff


  // opencl launch related stuff

  // Obtain the number of platforms available.
  err = clGetPlatformIDs(0, NULL, &num_platforms);
  if (err != CL_SUCCESS) process_error(err, "clGetPlatformIDs");

  platforms = malloc(sizeof(cl_platform_id) * num_platforms);

  // Obtain the list of platforms available.
  err = clGetPlatformIDs(num_platforms, platforms, NULL);
  if (err != CL_SUCCESS) process_error(err, "clGetPlatformIDs");

  // Get specific information about the OpenCL platforms.
  for (i = 0; i < num_platforms; i++) {

    // Obtain the list of devices available on the platform i.
    err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
    if (err != CL_SUCCESS) process_error(err, "clGetDeviceIDs");

    devices = malloc(sizeof(cl_device_id) * num_devices);
    err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);
    if (err != CL_SUCCESS) process_error(err, "clGetDeviceIDs");

  }

  // Create a context.
  context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateContext");
  
  // Create a command-queue on a specific device.
  // Here we are using the first device of the context.
  device = devices[0];
  command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateCommandQueue");

  // Create buffer objects. They will be the arguments to the kernel.
  flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
  d_A = clCreateBuffer(context, flags, sizeA_bytes, h_A, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateBuffer");

  d_B = clCreateBuffer(context, flags, sizeB_bytes, h_B, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateBuffer");

  d_nA = clCreateBuffer(context, flags, sizeof(size_t), &nA, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateBuffer");

  flags = CL_MEM_WRITE_ONLY;
  d_C = clCreateBuffer(context, flags, sizeC_bytes, h_C, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateBuffer");

  // Read program text file and place content into buffer (char *)
  program_filename = "matmult.cl";

  fd_program = fopen(program_filename, "r");
  if (fd_program == NULL) {
    perror("Error opening file");
    exit(EXIT_FAILURE);
  }

  fseek(fd_program, NO_OFFSET, SEEK_END);
  program_size = ftell(fd_program);
  rewind(fd_program);
  program_buffer = (char *) malloc(program_size + 1);
  fread(program_buffer, sizeof(char), program_size, fd_program);
  program_buffer[program_size] = '\0';
  fclose(fd_program);

  // Create a program object for the context, and load the source code specified by the text strings in the strings array into the program object.
  program = clCreateProgramWithSource(context, 1, (const char **) &program_buffer, &program_size, &err);
  if (err != CL_SUCCESS) process_error(err, "clCreateProgramWithSource");

  free(program_buffer);

  // Once readed the program, it is necessary to build it.
  // If device_list is NULL value, the program executable is built for all devices associated with program for which a source or binary has been loaded.
  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  if (err != CL_SUCCESS) {
    perror("clBuildProgram");
    char *build_log;
    // Find size of log and print to std output
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &param_value_size_ret);
    build_log = (char *) malloc(param_value_size_ret + 1);
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, param_value_size_ret, build_log, NULL);
    build_log[param_value_size_ret] = '\0';
    printf("CL_PROGRAM_BUILD_LOG: %s\n", build_log);
    free(build_log);
    exit(EXIT_FAILURE);
  }

  // Release the resources allocated by the OpenCL compiler. They are not needed any more.
  clUnloadCompiler();

  // Create a kernel object. (A kernel is a function declared in a program.)
  kernel = clCreateKernel(program, "matmult", &err); 
  if (err != CL_SUCCESS) process_error(err, "clCreateKernel");

  // Set the argument values for the specific arguments of the kernel.
  // Here we use the buffer objects created before.
  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_A);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_B);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_C);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_nA);
  if (err != CL_SUCCESS) process_error(err, "clSetKernelArg");


  // Specify the total number of threads to be launched on each dimension.
  global_work_size = (size_t *) malloc(DIMENSIONS * sizeof(size_t));

  global_work_size[0] = nA;
  global_work_size[1] = mB;

  // Enqueue the command to execute the kernel on the device.
  // As we do not set local_work_size, we let the implementation decide the amount of threads to run on each dimension of each work-group.
  work_dim = DIMENSIONS;
  err = clEnqueueNDRangeKernel(command_queue, kernel, work_dim, NULL, global_work_size, NULL, 0, NULL, &event);
  if (err != CL_SUCCESS) process_error(err, "clEnqueueNDRangeKernel");

  // Wait on the host thread for the command identified by the event object to complete.
  err = clWaitForEvents(1, &event);
  if (err != CL_SUCCESS) process_error(err, "clWaitForEvents");

  // Declaro las variables de tiempo
  cl_ulong time_start, time_end;
  double total_time;

  // Se extrae la información del tiempo
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
  total_time = time_end - time_start;
  printf("\nExecution time in milliseconds = %0.3f ms\n", (total_time / 1000000.0) );

  // Enqueue command to read asynchronously (CL_FALSE) from the buffer object to host memory.
  // We want to get the result computed by the kernel.
  err = clEnqueueReadBuffer(command_queue, d_C, CL_FALSE, NO_OFFSET, sizeC_bytes, h_C, 0, NULL, &event);
  if (err != CL_SUCCESS) process_error(err, "clEnqueueReadBuffer");

  // Wait on the host thread for the command identified by the event object to complete.
  err = clWaitForEvents(1, &event);
  if (err != CL_SUCCESS) process_error(err, "clWaitForEvents");

  // Esperar la finalizacion de todos los threads
  clFinish(command_queue);

  // Print something.
  //for (i = 0; i < sizeC ; i++) 
    //printf("%lf ", h_A[i]);
  //printf("\n");
  //for (i = 0; i < sizeC ; i++) 
    //printf("%lf ", h_B[i]);
  //printf("\n");
  //for (i = 0; i < sizeC ; i++) 
    //printf("%lf ", h_C[i]);

  // Deallocate resources
  clReleaseKernel(kernel);
  clReleaseMemObject(d_A);
  clReleaseMemObject(d_B);
  clReleaseMemObject(d_C);
  clReleaseMemObject(d_nA);
  clReleaseCommandQueue(command_queue);
  clReleaseProgram(program);
  clReleaseContext(context);

  // Free memory.
  free(global_work_size);
  free(h_A);
  free(h_B);
  free(h_C);
  free(devices);
  free(platforms);


  return EXIT_SUCCESS;
}
