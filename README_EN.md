This is an automatic translation, may be incorrect in some places. See sources and examples!

# upid
Easy Universal Library of PID-regulator with optional superstructures over the algorithm
- Limiting output
- Choosing the direction of regulation
- proportional to the error or input change
- differential error or entrance
- Integral component:
- restriction of the integral Back Calculation
- restriction of integration on the saturation of the exit
- The reset integrated when reached the setup

## compatibility
Compatible with all platforms

## Content
- [use] (#usage)
- [versions] (#varsions)
- [installation] (# Install)
- [bugs and feedback] (#fedback)

<a id = "USAGE"> </A>

## Usage
## upid
Standard class.The modes can be changed during the program, the coefficients can be written and read directly in the variables `kp`,` ki`, `kd`.

`` `CPP
upid (uint8_t cfg = 0, uint16_t dt = 30);

Float KP = 0, Ki = 0, KD = 0;
Float KBC = 0;
Float setpoint = 0;
Float Integral = 0;
Float outmax = 255;
Float outmin = 0;

Float getkp ();
Float getki ();
Float getkd ();

VOID Setkp (Float P);
VOID Setki (Float I);
VOID Setkd (Float D);

// Install config
VOID setconfig (Uint8_T NCFG);

// Turn on the flag
VOID setmode (upid_cfg mode);

// Turn off the flag
Void Clearmode (Upid_cfg Mode);

// Install the period of work in MS
VOID setdt (uint16_t ms);

// Calculate (call with a given period).Will return the way out
Float Compute (Float Input);
`` `

## upidfast
Fast and light version: the mode is set in the template at the compilation stage, and the coefficients can only be changed through SET-GET.This version is 20 μs faster - 67 versus 87 μs (AVR 16 MHZ, standard settings).

Also, when creating, you can specify a third-party implementation of `float`, for example, a number with a fixed point (library [fixed] (https://github.com/gyverlibs/fixed)), which will facilitate the code for 1-2 KB and accelerate to 40 μs (for processors without FPU).

`` `CPP
Upidfast <uint8_t cfg = 0, typename float_t = float> (uint16_t dt = 30);

Float_t KBC = 0;
Float_t setpoint = 0;
Float_t Integral = 0;
Float_t outmax = 255;
Float_t outmin = 0;

Float_t getkp ();
Float_t getki ();
Float_t getkd ();

VOID Setkp (float_t p);
VOID Setki (float_t i);
VOID Setkd (float_t d);

// Install the period of work in the MS (performed for a long time!)
VOID setdt (uint16_t ms);

// Calculate (call with a given period).Will return the way out
Float_t Compute (Float_T Input);
`` `

### Settings
The regime of the regulator is built from constants separated by `|`, for example:

`` `CPP
Upid pid (d_input | p_measure | pid_reverse);
pid.setconfig (d_input | p_measure);
pid.clearmode (p_measure);

Upidfast <d_input |P_measure |PID_REVERSE> PIDFAST;
`` `

## proportional
One option out of two:

- `p_error` (by default) - proportional to the error.The classic version of the P-component
- `p_measure` - proportional to the entrance (depends on` d_input`) or error (depending on `d_error`).Changes the logic of the regulator, well suited for integrating processes - for example, the regulator controls the speed (signal to the motor), and the input is given a position (encoder on this motor)

### Integral
You can choose in any combination, but `` I_BACK_CALC` is analCranberries `i_saturate` and turning them together does not make sense:

- `I_SATURATE` - *Conditional integration *, disconnecting integration when saturating the output (` Outmax` and `Outmin` must be configured)
- `I_Back_calc` - *Back Calculation *, a smart restriction of integration when saturating the output (` outmax` and `Outmin` should be configured).It works slower than `` _SATURATE`.The intensity is regulated by the CoEF -M `KBC` - first you need to choose it equal to` ki`, then change and look behind the reaction of the system
- `I_Reset` - automatic discharge integral when reached the setup

> By "saturation of the exit", it means his departure from the `` (out .. outmax) `` `` `Outmar, in ordinary PID, the integral sum begins to grow uncontrolled in one direction and then returns for a long time - the so -called*Windup*.The library offers several options for solving this problem

### differential
One option out of two:

- `D_ERROR` (by default) - differentiation of error.The classic version of D-component
- `d_input` - the differentiation of the entrance.The less impact of noise, the complete absence of casting in D when changing the setting.If the setting often changes and you need more dynamic behavior of the system in response to this - it is better to use `d_error`

### Direction
One option out of two:

- `pid_forward` (default) - direct regulator, the regulator’s output should increase the entrance (for example, heating power control, the temperature to the regulator input)
- `pid_reverse` - the return regulator, the regulator output should reduce the entrance (for example, cooling power control, the temperature to the regulator input)

### Examples
`` `CPP
#include <upid.h>

Consta int dt = 30;
Upid pid (d_input | _saturate);
// upidfast <d_input |I_SATURATE> pid;

VOID setup () {
// pid.kp = 10;
// pid.ki = 20;
// pid.kd = 5;

Pid.Setkp (10);
Pid.Setki (20);
Pid.Setkd (5);

// pid.kbc = 0.1;
Pid.Setdt (DT);
Pid.outmax = 255;
pid.outmin = -255;

Pid.Setpoint = setting;
}

VOID loop () {
Float Result = Pid.compute
apply (result);
DELAY (DT);
}
`` `

<a ID = "Versions"> </a>

## versions
- V1.0

<a id = "Install"> </a>
## Installation
- The library can be found by the name ** upid ** and installed through the library manager in:
- Arduino ide
- Arduino ide v2
- Platformio
- [download library] (https://github.com/gyverlibs/upid/archive/refs/heads/main.zip). Zip archive for manual installation:
- unpack and put in * C: \ Program Files (X86) \ Arduino \ Libraries * (Windows X64)
- unpack and put in * C: \ Program Files \ Arduino \ Libraries * (Windows X32)
- unpack and put in *documents/arduino/libraries/ *
- (Arduino id) Automatic installation from. Zip: * sketch/connect the library/add .Zip library ... * and specify downloaded archive
- Read more detailed instructions for installing libraries[here] (https://alexgyver.ru/arduino-first/#%D0%A3%D1%81%D1%82%D0%B0%D0%BD%D0%BE%D0%B2%D0%BA%D0%B0_%D0%B1%D0%B8%D0%B1%D0%B8%D0%BE%D1%82%D0%B5%D0%BA)
### Update
- I recommend always updating the library: errors and bugs are corrected in the new versions, as well as optimization and new features are added
- through the IDE library manager: find the library how to install and click "update"
- Manually: ** remove the folder with the old version **, and then put a new one in its place.“Replacement” cannot be done: sometimes in new versions, files that remain when replacing are deleted and can lead to errors!

<a id = "Feedback"> </a>

## bugs and feedback
Create ** Issue ** when you find the bugs, and better immediately write to the mail [alex@alexgyver.ru] (mailto: alex@alexgyver.ru)
The library is open for refinement and your ** pull Request ** 'ow!

When reporting about bugs or incorrect work of the library, it is necessary to indicate:
- The version of the library
- What is MK used
- SDK version (for ESP)
- version of Arduino ide
- whether the built -in examples work correctly, in which the functions and designs are used, leading to a bug in your code
- what code has been loaded, what work was expected from it and toCranberries he works in reality
- Ideally, attach the minimum code in which the bug is observed.Not a canvas of a thousand lines, but a minimum code