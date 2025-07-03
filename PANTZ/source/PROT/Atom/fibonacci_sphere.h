/* Created by the clay at Auburn University
 *
 * This file is intended to be included by the Atom.h header file. It contains
 * the fibbonaci sphere building method of the Atom class. */

// Confirm that the Atom class has been declared and is actively being loaded
#ifndef Atom_Loading_Status
#error Atom methods must be included from the Atom.h header file
#endif 

// the sphere around the atom (for interaction calculations)
vector<vector<float>> PROT::Atom::fibonacci_sphere (size_t N, float probe_radius){
    // the vector of points that represent the sphere
    vector<vector<float>> points;
    // the radius of this sphere
    float R = probe_radius + this->lj_sigma();
    // the golden angle (angle found from golden ratio cirumference of a circle)
    float golden_ratio = (1 + sqrt(5)) / 2;
    // loop through the number of points
    for (size_t k = 0; k < N; k++) {
        float theta = 2 * M_PI * k / golden_ratio;
        float phi = acos(1 - 2 * (k + 0.5) / N);
        // the x and z coordinates
        float x = sin(phi) * cos(theta);
        float y = sin(phi) * sin(theta);
        float z = cos(phi);
        // scale the coordinates by the radius and atom location
        points.push_back({(x*R + m_coors[0]), (y*R + m_coors[1]), (z*R + m_coors[2])});
    }
    // Return the sphere
    return points;
}

// a function to get the radius of the atom (lennard jones sigma)
float PROT::Atom::lj_sigma() const {
    // Radii that Varun used
    // https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html
    if (m_element == "H") {
        return 1;
    }
    else if (m_element == "C") {
        return 1.7;
    }
    else if (m_element == "N") {
        return 1.625;
    }
    else if (m_element == "O") {
        return 1.5;
    }
    else if (m_element == "S") {
        return 1.782;
    } else {
        string error = "Atom name not found in sigma map: " + m_element;
        cout<<error<<endl;
        throw error;
    }
}

