/* Created by the Pantazes Lab at Auburn University. 
 *
 * This file contains the definition of the PDB class and includes the header
 * files that implement its methods. */

// Use a header guard to prevent this file from being included in a compiled
// program multiple times
#ifndef Proteins_PDB_Guard
#define Proteins_PDB_Guard 1

// Make sure this file is being included from the Proteins.h header file
#ifndef Proteins_Loading_Status
#error PDB.h must be included by Proteins.h
#endif

// Confirm that the Structure class has been loaded
#ifndef Proteins_Structure_Guard
#error PDB.h must be included after Structure.h
#endif

// The PDB class organizes the information from a PDB file into a useable form
class PROT::PDB {

    // The information stored in the class is private
    private:
        // The name of the pdb file
        string m_name;
        // The folder the file is located in
        string m_folder;
        // The contents of that file
        vector<string> m_lines;
        // The type of experiment used to generate the information
        string m_type;
        // The resolution of the experimental data
        double m_resolution;
        // The Proteins in the file
        vector<Protein> m_proteins;
        // The structures they group into
        vector<Structure> m_structures;
        // Whether or not the PDB file is obsolete and shouldn't be used
        bool m_obsolete;
        // Whether or not the PDB file is a theoretical structure
        bool m_theoretical;

    // Private functions that control behaviour of the class
    private:
        // Copy information from another instance of this class
        void copy (const PDB *);
        // Load the contents of the file
        void load ();
        // Identify the experiment type
        void identify_type();
        // Identify the experiment's resolution
        void identify_resolution();
        // Identify residues that are missing using the REMARK 465 lines
        void identify_missing_residues (vector<vector<Residue> >&);
        // identify Atoms that are missing using the REMARK 470 lines
        void identify_missing_atoms ();
        // Check to see whether or not the PDB file is obsolete or a
        // theoretical model. 
        void check_file_status ();
        // Use the SEQRES lines to initialize vectors of Residues that will
        // eventually be turned into Proteins
        void initialize_Residues (vector<vector<Residue> >&);
        // Integrate missing residues and make Proteins
        void integrate_missing_residues (Protein&, vector<Residue>&);
        // Collect the Atoms that provide structural details about the
        // Proteins
        void collect_Atoms (vector<vector<Atom> >&);
        // Validate the SEQRES information
        void validate_seqres (Protein *, vector<Residue>&);
        // Construct the proteins
        void construct_Proteins ();
        // Create the Structures
        void create_Structures ();

    // The public interface of the PDB class
    public:
        // The class constructor
        PDB (const string&, const string);
        // Copy construction and assignment
        PDB (const PDB& other) {copy(&other);}
        PDB (const PDB * other) {copy(other);}
        void operator= (const PDB& other) {copy(&other);}
        void operator= (const PDB * other) {copy(other);}
        // Access to the class information
        string name () const {return m_name;}
        string folder () const {return m_folder;}
        size_t lines () const {return m_lines.size();}
        string line (const size_t) const;
        string type () const {return m_type;}
        double resolution () const {return m_resolution;}
        size_t proteins () const {return m_proteins.size();}
        Protein * protein (const size_t);
        size_t structures () const {return m_structures.size();}
        Structure * structure (const size_t);
        bool obsolete () const {return m_obsolete;}
        bool theoretical () const {return m_theoretical;}
        // A string representation of the PDB file's information
        string str () const;
        // set the elements of the atoms in the PDB file
        void set_elements ();

    // End the class definition
};

// Define a preprocessor variable so the class methods can confirm they are
// being loaded when appropriate
#define PDB_Loading_Status 1

// Include the PDB header files
#include "PDB/copy.h"
#include "PDB/load.h"
#include "PDB/identify.h"
#include "PDB/file_status.h"
#include "PDB/residues.h"
#include "PDB/collect_atoms.h"
#include "PDB/validate.h"
#include "PDB/create.h"
#include "PDB/constructor.h"
#include "PDB/line.h"
#include "PDB/protein.h"
#include "PDB/structure.h"
#include "PDB/str.h"
#include "PDB/set_elements.h"

// Undefine the preprocessor variable
#undef PDB_Loading_Status

// End the if statement from the start of the file
#endif
