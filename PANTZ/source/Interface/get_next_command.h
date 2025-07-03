/* Created by the Pantazes Lab at Auburn University.
 * 
 * This file gets the next command line from the provided instructions. */

// This file is intended to be included directly from the Interface.h header
// file - check that that is the case
#ifndef PANTZ_Interface_Loading_Status
#error Interface methods must be included by Interface.h
#endif

// Implement the method
size_t Interface::get_next_command (vector<string>& parts, 
                                    size_t Index,
                                    const bool throwError = true) {
    // Repeat this until either the end of the commands are reached or the next
    // non-comment / blank line is found
    while (Index < m_all_contents.size()) {
        // Get the line with out any comments or whitespace
        string line = remove_comments(m_all_contents[Index]);
        // If the line is empty, increment the index and continue the loop
        if (line.size() == 0) {Index += 1; continue;}
        // Since the line isn't empty, use it to fill in the parts. This
        // function can throw an error. If it does, that is because of a
        // formatting issue in the identified commmand. That error SHOULD always
        // be thrown by this function.
        create_parts (parts, line, Index);
        // Break the search
        break;}
    // If the end of the file has been reached and that wasn't supposed to
    // happen, throw an error
    if ((Index >= m_all_contents.size()) && (throwError)) {
        string error = "The end of the instruction file was reached "
                       "unexpectedly.\n";
        throw PANTZ_error (error);}
    // Return the index of the command that was checked most recently
    return Index;
}
