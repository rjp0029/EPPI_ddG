/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the declaration and implementation of an error class that
 * should be used whenever errors occur in the PANTZ codes. */

// Use a header guard to make sure the file is only included in a compiled
// program a single time
#ifndef PANTZ_Error_Guard
#define PANTZ_Error_Guard 1

// Include the Text.h header file, which includes the string header. It also
// declares that the standard namespace is being used
#include "Text.h"
// Include the exception header file
#include <exception>

// Declare the class
class PANTZ_error : public exception {

    // The information in the class that is private
    private:
        // The error message
        string m_message;

    // The public interface of the class
    public:
        // The cluster compiler threw an error about how the PANTZ_error
        // destructor (which was not defined) was overwriting the exception
        // destructor with a looser throw specifier. This code has been
        // introduced for the sole purpose of correcting that compile error
        ~PANTZ_error () throw () {m_message = "";}

        // Initialize from a string
        PANTZ_error (const string& input) {m_message = input;}

        // Initialize from a character pointer. Given the programming habits of
        // Dr. Pantazes, this option is probably never used.
        PANTZ_error (const char * input) {m_message = input;}

        // Initialize from another instance of the class and a string 
        PANTZ_error (const PANTZ_error& other, const string& input) {
            // Start the message out as an empty string
            m_message = "";
            // If the input string isn't empty, start the message with that
            if (input.size() > 0) {
                m_message += input;
                // If the message isn't ending with an end line character, add
                // one
                if (!Text::endline(m_message)) {m_message += "\n";}}
            // Add the previous error's message
            m_message += other.m_message;}

        // A default initialization
        PANTZ_error () {
            m_message = "An error occurred in the PANTZ calculations.\n";}

        // The typical 'what' method of an exception
        virtual const char* what () const throw() {
            const char* output = m_message.c_str();
            return output;}

// End the class definition
};

// End the header guard from the start of the file
#endif
