#ifndef __CALLBACK_H__
#define __CALLBACK_H__

template <typename T, typename ... S>
struct Callback
{
    virtual T operator()(S... args) = 0;
    virtual ~ Callback () { }
};

#endif
