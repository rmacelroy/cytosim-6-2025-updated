// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// This file create by FJN on 27/9/2008.

#ifndef TIME_DATE_H
#define TIME_DATE_H

#include <ctime>

/// A set of functions related to time
/**
 Functions to get (real) time and date, from the C-standard library.
 */
namespace TimeDate
{
    /// set current date in short format, `buf` should be 26 character long or more
    void get_date(char * buf, size_t buf_size);
    
    /// set current date in short format, `buf` should be 26 character long or more
    void get_date(char * buf, size_t buf_size, bool no_year);

    /// using a local char[] buffer to call get_date()
    char const* date_string();

    
    /// approximately the number of days since Jan 1 2000
    int days_since_2000();
    
    /// number of seconds since Jan 1 1970, 0h00
    time_t seconds_since_1970();
    
    /// number of seconds since Jan 1 2000, 0h00
    time_t seconds_since_2000();

    /// year
    int year();

    /// day of the year (0-365)
    int day_of_the_year();
    
    /// hour of the day (0-23)
    int hours_today();

    /// number of seconds since midnight
    double seconds_today();
    
    /// number of milliseconds since midnight
    double milliseconds();
 
    /// number of microseconds
    double microseconds();
    
    /// return CPU time in seconds and update argument to current time
    double processor_time(clock_t&);

}

#endif
