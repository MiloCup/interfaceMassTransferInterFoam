//D =  D_Int
//volScalarField D ;

// set internal D field
forAll(D, cellI)
   {

   if (alpha1[cellI] < scalar(0.95) &&   alpha1[cellI] > scalar(0.05))
   {
      D[cellI] = D_Int.value();

   }
   else
      if
      (
        alpha1[cellI] > scalar(0.95)
      )

         {
          D[cellI] = D_Water.value();
         }
      else
//      (
//        alpha1[cellI] < scalar(0.05)
//      )

         {
          D[cellI] = D_Air.value();
         }
   }

// set boundary field
forAll(D.boundaryField(), patchI)
      {
         // Alpha1 is water.  This checks that alpha1 is water enough and sets
         forAll(D.boundaryField()[patchI], faceI)
         {

            if (alpha1.boundaryField()[patchI][faceI]  < scalar(0.95) &&   alpha1.boundaryField()[patchI][faceI]  > scalar(0.05))
            {
               D.boundaryFieldRef()[patchI][faceI] = D_Int.value();
            }
            else
               if
                  (
                     alpha1.boundaryField()[patchI][faceI] > scalar(0.95)

                  )
               {
                D.boundaryFieldRef()[patchI][faceI] = D_Water.value();
               }
               else
//                  (
//                     (alpha1.boundaryField()[patchI][faceI]) < scalar(0.05)
//                  )
               {
                D.boundaryFieldRef()[patchI][faceI] = D_Air.value();
               }
         }
      }
/// update coupled boundaries
D.correctBoundaryConditions();

//record updated D
if(runTime.outputTime())
{
D.write();

}


fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(D, T)
             ==
                fvOptions(T)
            );

TEqn.relax();
fvOptions.constrain(TEqn);
TEqn.solve();
fvOptions.correct(T);
