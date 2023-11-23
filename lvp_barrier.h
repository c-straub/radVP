#ifndef BARRIER_H_INCLUDED
#define BARRIER_H_INCLUDED


#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;

int number_of_threads_const = 200; //Number of threads you wish to use (maximal usable number probably depends on your computer/hardware, but can be larger than the number of cores)



mutex muter;



/////////////////////////////////////////////////
/// Helper for parallelization and barriers
//  The following is needed to implement parallelization with the std thread package in C++11, in particular barriers (which make the threads wait for each other).
//  By Christopher Straub (with help from stackoverflow), 11.08.2019
/////////////////////////////////////////////////

class Barrier //Helper class for barriers
{
    private:
        mutex _mutex;
        condition_variable _cv;
        size_t _count;
    public:
        explicit Barrier(size_t count) : _count(count) { }
        void wait()
        {
            unique_lock<mutex> lock(_mutex);
            if (--_count == 0)
                _cv.notify_all();
            else
                _cv.wait(lock, [this] { return _count == 0; });
        }
};
vector<Barrier*> barriers; //Global array containing every used barrier object (as a pointer)
int counter_of_used_barriers = 0;
int counter_of_function_calls = 0;

inline void wait_for_all_threads(int barrier_counter)
{
    if(barrier_counter < barriers.size())
        barriers[barrier_counter]->wait();
    else
    {
        {
            lock_guard<mutex> lock(muter); //Makes sure that only one thread creates the new Barrier at a time;
            while(!(barrier_counter < barriers.size()))
                barriers.push_back(new Barrier(number_of_threads_const)); //Note that we need the number of threads here, therefore it was fixed at the very top of this file
        } //NOTE that the mutex object will be unlocked if the deconstructor of the lock_guard object is called; which is here, since we leave the scope where it has been defined...
        barriers[barrier_counter]->wait();
    }
}

inline void wait_for_all_threads() //Calls the function above and has an internal counter
{
    {
        lock_guard<mutex> lock(muter);
        counter_of_function_calls++;
        if(counter_of_function_calls > (counter_of_used_barriers+1)*number_of_threads_const)
        {
            counter_of_used_barriers++;
            if(counter_of_used_barriers > 10000) //prevents the barrier array from getting infinitely large... (since "counter_of_function_calls" overflowed in the past...)
            {
                barriers.clear();
                counter_of_used_barriers = 0;
                counter_of_function_calls = 1;
            }
        }
    }
    wait_for_all_threads(counter_of_used_barriers);
}

///IMPORTANT NOTE: For the functions above to work properly, one has to ALWAYS call the wait_for_all_threads function with every thread (i.e., always with number_of_threads_const threads)





#endif // BARRIER_H_INCLUDED
