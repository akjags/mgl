#ifdef documentation
=========================================================================

       program: writeDigPort.c
            by: justin gardner
          date: 06/26/09
  adapted from: NI-DAQmx Base example file writeDigPort.c -- made into
                a mex file and set to only init task once. Note that
                if you need to write from multiple ports, you may want to 
                modify this file to be more efficient -- right now it 
                will have to close and open each port to read from different
                ports which will make it slow. If you are reading from a 
                single port then you should be ok.


=========================================================================
#endif
/*********************************************************************
*
* C Example program:
*    readDigPort.c
*
* Example Category:
*    DI
*
* Description:
*    This example demonstrates how to read values from a digital input
*    channel.
*
* Instructions for Running:
*    1. Select the digital ports on the DAQ device to be read.
*
* Steps:
*    1. Create a task.
*    2. Create a Digital Input channel.
*    3. Call the Start function to start the task.
*    4. Read the digital uInt8 array data. This read function
*       reads a single sample of digital data on demand, so no timeout
*       is necessary.
*    5. Call the Clear Task function to clear the Task.
*    6. Display an error if any..
*
*********************************************************************/
/////////////////////////
//   include section   //
/////////////////////////
#include "/Applications/National Instruments/NI-DAQmx Base/includes/NIDAQmxBase.h"
#include <stdio.h>
#include "mex.h"

////////////////////////
//   define section   //
////////////////////////
#define DAQmxErrChk(functionCall) { if( DAQmxFailed(error=(functionCall)) ) { goto Error; } }

//////////////////////
//   function decls //
//////////////////////
int startTask(int portNum);
void stopTask(int portNum);

//////////////////////////
//   static variables   //
//////////////////////////
static TaskHandle  taskHandle = 0;
static int portNum = 1;

///////////////////////
//   main function   //
///////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Task parameters
   int32       error = 0;
   char        errBuff[2048];

   // Channel parameters
   char chan[256];

   // write variables
   int32       written;

   // pointer to output data
   double *outptr;
   uInt32 val, channelNum = 2;

   // set up return value
   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   outptr = mxGetPr(plhs[0]);

   // check input arguments
   if (nrhs < 1) {
     mxArray *callInput[] = {mxCreateString("writeDigPort")};
     mexCallMATLAB(0,NULL,1,callInput,"help");
     *outptr = (double)0;
     return;
   }

   // get the channel number
   if (nrhs == 2) {
     if (mxGetPr(prhs[1]) != NULL) {
       // get the input port number
       int inputPortNum = (int) *mxGetPr( prhs[1] );
       // see if we are being asked to shutdown
       if (inputPortNum == -1) {
	 stopTask(portNum);
	 return;
       }
       // portNum == -2 means to print info
       else if (inputPortNum == -2) {
	 if (taskHandle == 0)
	   mexPrintf("(writeDigPort) No write task is currently initialized\n");
	 else
	   mexPrintf("(writeDigPort) Write task open on port %i\n",portNum);
	 return;
       }
       // have to restart task if we are reading from a different port
       if (inputPortNum != portNum) stopTask(portNum);
       portNum = inputPortNum;
     }
   }

   // Create Digital Output (DO) Task and Channel, only if the task is not initialized
   if (taskHandle == 0) {
     if (startTask(portNum) == 0) {
       // if we were unsuccessful then return 
       *outptr = (double)0;
       return;
     }
   }

   // get input value
   val = (uInt32)*mxGetPr(prhs[0]);
   // DAQmxBaseWriteDigitalU8 
   DAQmxErrChk (DAQmxBaseWriteDigitalU32(taskHandle,1,1,10.0,DAQmx_Val_GroupByChannel,&val,&written,NULL));

   // return success
   *outptr = (double)1;

   return;

 Error:

   if (DAQmxFailed (error))
     DAQmxBaseGetExtendedErrorInfo (errBuff, 2048);

   if (taskHandle != 0)
   {
      DAQmxBaseStopTask (taskHandle);
      DAQmxBaseClearTask (taskHandle);
   }

   if ((error) && (error != -200220))
     mexPrintf ("(writeDigPort) DAQmxBase Error %d: %s\n", error, errBuff);

   *outptr = (double)-1;

   return;
}

///////////////////
//   startTask   //
///////////////////
int startTask(int portNum)
{
  // Error variables
  int32       error = 0;
  char        errBuff[2048];
  
  // Setup the channel parameter
  char  chan[] = "Dev1/port1";
  if (portNum == 0)
    chan[9] = '0';
  else if (portNum == 2)
    chan[9] = '2';
  else if (portNum != 1) {
    mexPrintf("(writeDigPort) Unrecogonized port number %i. Should be 0-2\n",portNum);
    return 0;
  }
		   
  // Create the task
  DAQmxErrChk (DAQmxBaseCreateTask ("", &taskHandle));
  DAQmxErrChk (DAQmxBaseCreateDOChan(taskHandle,chan,"",DAQmx_Val_ChanForAllLines));
   
  // Start Task (configure port)
  DAQmxErrChk (DAQmxBaseStartTask (taskHandle));

  // return success
  return 1;

 Error:

   if (DAQmxFailed (error))
     DAQmxBaseGetExtendedErrorInfo (errBuff, 2048);

   if (taskHandle != 0) stopTask(portNum);

   // output error, but only if it is not device idnetifier is invalid
   // since this happens when you simply don't have a card in the
   // computer
   if (error)
     if (error != -200220)
       mexPrintf ("(writeDigPort) DAQmxBase Error %d: %s\n", error, errBuff);
     else
       mexPrintf ("(writeDigPort) DAQmxBase Error %d: %s\n", error, errBuff);
       
   return 0;
}

//////////////////
//   stopTask   //
//////////////////
void stopTask(int portNum)
{
  if (taskHandle != 0) {
    // stop task
    DAQmxBaseStopTask (taskHandle);
    DAQmxBaseClearTask(taskHandle);
  }
  taskHandle = 0;
}
