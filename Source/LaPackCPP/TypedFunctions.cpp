#pragma once

/** For some of the functions in LibF2C, Blas and LaPack, it didn't make sense to templatize them.
These are all put into this file because it may be problematic in certain contexts to have mixed 
templated and typed code in a single implementation file - if you include the same implementation
file in two different unity-build compilation units, you will get "multiple definition" linker 
errors. You may have to include the same implementation file into different compilation unit, if
these compilation units instantiate different templates or the same template for different 
datatypes. I know - it's a mess :-( */

namespace LaPackCPP {

//=================================================================================================
// LibF2C



//=================================================================================================
// Blas


//=================================================================================================
// XBlas


//=================================================================================================
// LaPack







}