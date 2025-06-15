#pragma once
#include <stdint.h>

// uPID_cfg
enum uPID_cfg : uint8_t {
    P_ERROR = 0,             // пропорционально ошибке
    P_MEASURE = (1 << 0),    // пропорционально вводу D_INPUT или ошибке D_ERROR. Для интегрирующих процессов
    I_BACK_CALC = (1 << 1),  // Back Calculation, ограничение интегрирования регулируется коэф-м Kbc (по умолч. выбирать равным Ki)
    I_SATURATE = (1 << 2),   // Conditional Integration, отключение интегрирования при насыщении выхода
    I_RESET = (1 << 3),      // сброс интеграла при достижении уставки
    D_INPUT = (1 << 4),      // дифференциально входу
    D_ERROR = 0,             // дифференциально ошибке
    PID_FORWARD = 0,         // прямое направление
    PID_REVERSE = (1 << 5),  // обратное направление
};

// uPID
class uPID {
   public:
    uPID(uint8_t cfg = 0, uint16_t dt = 30) : cfg(cfg) {
        setDt(dt);
    }

    float Kp = 0, Ki = 0, Kd = 0;
    float Kbc = 0;
    float setpoint = 0;
    float integral = 0;
    float outMax = 255;
    float outMin = 0;

    inline float getKp() { return Kp; }
    inline float getKi() { return Ki; }
    inline float getKd() { return Kd; }

    inline void setKp(float p) { Kp = p; }
    inline void setKi(float i) { Ki = i; }
    inline void setKd(float d) { Kd = d; }

    // установить конфиг
    inline void setConfig(uint8_t ncfg) { cfg = ncfg; }

    // включить флаг
    inline void setMode(uPID_cfg mode) { cfg |= mode; }

    // выключить флаг
    inline void clearMode(uPID_cfg mode) { cfg &= ~mode; }

    // установить период работы в мс
    void setDt(uint16_t ms) {
        _dt = ms / 1000.0f;
        _dt_i = 1.0f / _dt;
    }

    // вычислить (вызывать с заданным периодом). Вернёт выход
    float compute(float input) {
        float error = setpoint - input;
        float deriv = (cfg & D_INPUT) ? (_prev - input) : (error - _prev);
        _prev = (cfg & D_INPUT) ? input : error;

        if (_first) _first = false, deriv = 0;

        if (cfg & PID_REVERSE) {
            error = -error;
            deriv = -deriv;
        }

        if ((cfg & I_RESET) && ((!(cfg & PID_REVERSE) && input >= setpoint) || ((cfg & PID_REVERSE) && input <= setpoint))) integral = 0;
        if (cfg & P_MEASURE) integral += Kp * deriv;

        float output = ((cfg & P_MEASURE) ? 0 : Kp * error) + integral + Kd * deriv * _dt_i;

        if (output >= outMax) {
            if ((cfg & I_BACK_CALC) && Kbc) error += (outMax - output) * Kbc;
            output = outMax;
            if ((cfg & I_SATURATE) && error > 0) return output;
        } else if (output <= outMin) {
            if ((cfg & I_BACK_CALC) && Kbc) error += (outMin - output) * Kbc;
            output = outMin;
            if ((cfg & I_SATURATE) && error < 0) return output;
        }

        integral += Ki * error * _dt;

        return output;
    }

   private:
    float _prev = 0;
    float _dt, _dt_i;
    bool _first = true;
    uint8_t cfg = 0;
};

// uPIDfast
template <uint8_t cfg = 0, typename float_t = float>
class uPIDfast {
    static constexpr bool _backCalc = cfg & I_BACK_CALC;
    static constexpr bool _satur = cfg & I_SATURATE;
    static constexpr bool _pMeas = cfg & P_MEASURE;
    static constexpr bool _rev = cfg & PID_REVERSE;
    static constexpr bool _dInput = cfg & D_INPUT;
    static constexpr bool _rst = cfg & I_RESET;

   public:
    uPIDfast(uint16_t dt = 30) {
        setDt(dt);
    }

    float_t Kbc = 0;
    float_t setpoint = 0;
    float_t integral = 0;
    float_t outMax = 255;
    float_t outMin = 0;

    inline float_t getKp() { return Kp; }
    inline float_t getKi() { return Ki / _dt; }
    inline float_t getKd() { return Kd * _dt; }

    inline void setKp(float_t p) { Kp = p; }
    inline void setKi(float_t i) { Ki = i * _dt; }
    inline void setKd(float_t d) { Kd = d / _dt; }

    // установить период работы в мс
    void setDt(uint16_t ms) {
        float_t r = (ms / 1000.0f) / _dt;
        Ki *= r;
        Kd /= r;
        _dt = ms / 1000.0f;
    }

    // вычислить (вызывать с заданным периодом). Вернёт выход
    float_t compute(float_t input) {
        float_t error = setpoint - input;
        float_t deriv = _dInput ? (_prev - input) : (error - _prev);
        _prev = _dInput ? input : error;

        if (_first) _first = false, deriv = 0;

        if (_rev) {
            error = -error;
            deriv = -deriv;
        }

        if (_rst && ((!_rev && input >= setpoint) || (_rev && input <= setpoint))) integral = 0;
        if (_pMeas) integral += Kp * deriv;

        float_t output = (_pMeas ? 0 : Kp * error) + integral + Kd * deriv;

        if (output >= outMax) {
            if (_backCalc && Kbc) error += (outMax - output) * Kbc;
            output = outMax;
            if (_satur && error > 0) return output;
        } else if (output <= outMin) {
            if (_backCalc && Kbc) error += (outMin - output) * Kbc;
            output = outMin;
            if (_satur && error < 0) return output;
        }

        integral += Ki * error;

        return output;
    }

   private:
    float_t Kp = 0, Ki = 0, Kd = 0;
    float_t _prev = 0;
    float_t _dt;
    bool _first = true;
};