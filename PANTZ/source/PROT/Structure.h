/* Created by the Pantazes Lab at Auburn University.
 *
 * This file contains the definitions of the Structure class, which is a
 * container of repeated Proteins in a PDB file. In a compiled program, it is
 * responsible for including the necessary header files for the class methods.
 * */

// Use a header guard to make sure this is only included once
#ifndef Proteins_Structure_Guard
#define Proteins_Structure_Guard 1

// Make sure the file is being included from Proteins.h
#ifndef Proteins_Loading_Status
#error Structure.h must be included by Proteins.h
#endif

// Confirm that the Protein class has been loaded
#ifndef Proteins_Protein_Guard
#error Structure.h must be included after Protein.h
#endif

// Define the class
class PROT::Structure {

    // The PDB class is a friend so it can access internal information
    friend class PROT::PDB;

    // The information stored in the class
    private:
        // Pointers to the Proteins in the Structure. Note that this class is
        // not responsible for storing the actual objects
        vector<Protein *> m_proteins;
        // The names of the structure
        vector<string> m_names;

    // Private methods for copying class content
    private:
        void copy (const Structure *);
        void copy (const Structure *, vector<Protein>&);

    // The public interface of the class
    public:
        // The default constructor
        Structure () {m_proteins.clear(); m_names.clear();}
        // Copy construction and assignment
        Structure (const Structure& other) {copy(&other);}
        Structure (const Structure * other) {copy(other);}
        void operator= (const Structure& other) {copy(&other);}
        void operator= (const Structure * other) {copy(other);}
        // Copy construction for use in PDB file copying
        Structure (const Structure& other, vector<Protein>& prots) {
            copy(&other, prots);}
        // The number of proteins in the Structure
        size_t proteins () const {return m_proteins.size();}
        // Access to a pointer
        PROT::Protein * protein (const size_t);
        // Naming information
        size_t names () const {return m_names.size();}
        string name (const size_t) const;
        // Store information
        void store_protein (Protein * ptr) {m_proteins.push_back(ptr);}
        void store_name (const string& text) {m_names.push_back(text);}
        // A string summary of the class information
        string str() const;

    // End the class definition
};

// Use a loading status to make sure headers are all included here and only here
#define Structure_Loading_Status 1
#include "Structure/copy.h"
#include "Structure/get.h"
#include "Structure/str.h"

// Undefine the structure loading status
#undef Structure_Loading_Status

// End the header guard from the start of the file
#endif
