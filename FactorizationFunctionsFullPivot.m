mA=[];

## Helper Function 
function isValid = validateEquation(mA,mb)
  isValid = false;
  if (rows(mA) != rows(mb)) 
    disp("A and b must have the same number of rows")
    return
  endif
  isValid = true;
endfunction
 

## Helper Function
#   Find the optimal pivot given permutation matrix P
#   m = P*A -- using indirect references 
function P = findPivot(m,P,n,pivot) 
  row = pivot;
  for i = pivot:n
    if (abs(m(P(row),pivot)) < abs(m(P(i),pivot))) 
      row = i;
    endif
  endfor
  
  if (row != pivot) 
    # Need to swap rows 
    temp = P(pivot);
    P(pivot) = P(row);
    P(row) = temp;
  endif
endfunction 


## Lab 1 - Function #1
function x = guassianElimination(A,b) 
  if (rows(b) == 1 && columns(b) != 1) 
    b = transpose(b);
  endif
  
  # Validation 
  if (!validateEquation(A,b))
    x = 0;
    return;
  endif
  
  m = [A b];
  r = rows(m);
  c = columns(m);
  
  # Setup the permutation tracking 
  P = 1:r;
  lowTri = eye(r);
  
  for pivot = 1:r
    # FULL PIVOTING 
    P = findPivot(m,P,r,pivot);
    
    # Create zeros in the rows other than the pivot 
    for row = 1:r 
      if (row == pivot) 
        continue; 
      endif
      rowFactor = (m(P(row),pivot) / m(P(pivot),pivot));
      m(P(row),:) = m(P(row),:) - (rowFactor*m(P(pivot),:));
    endfor
  endfor
  
  # The diagonal is not yet all ones - and the permutation needs to be accounted for
  for i = 1:r
    x(i) = m(P(i),c) / m(P(i),i);
  endfor
endfunction

## Lab 1 - Function 2
function [lowTri,upTri,P] = factorizationLR(A) 
  m = A;  
  r = rows(m);
  c = columns(m);
  
  # Setup the permutation tracking
  P = 1:r;
  lowTri = eye(r);
  
  for pivot = 1:r
    # FULL PIVOTING 
    P= findPivot(m,P,r,pivot);
    
    # Create zeros in the rows "below" the pivot 
    for row = (pivot+1):r
      lowTri(row,pivot) = (m(P(row),pivot) / m(P(pivot),pivot));
      m(P(row),:) = m(P(row),:) - (lowTri(row,pivot)*m(P(pivot),:));
    endfor
  endfor

  # Rearrange the rows so that they are in order based on the permutation matrix 
  upTri(:,:)=m(P,:);
endfunction 


## Lab 1 - Function 3 (part 1)
function c = forwardSubstitution(L,b) 
  if (rows(b) == 1 && columns(b) != 1) 
    b = transpose(b);
  endif
  
  # Validate that the rows are the same 
  if (!validateEquation(L,b))
    return; 
  endif
  
  r = rows(L);
  
  c(1) = b(1)/L(1,1);
  for i=2:r
    #otherTerms = 0;
    #for j=1:i-1
    #  otherTerms = otherTerms + (L(i,j)*c(j));
    #endfor
    #c(i) = (1/L(i,i))*(b(i) - otherTerms);
    c(i) = (1/L(i,i))*(b(i) - dot(L(i,1:i-1),c(1:i-1)));
  endfor
endfunction

## Lab 1 - Function 3 (part 2)
function x = backwardSubstitution(U,c) 
  r = rows(U);
  
  x(r) = c(r)/U(r,r);
  for i=(r-1):-1:1
    #otherTerms = 0;
    #for j=i+1:r
    #  otherTerms = otherTerms + (U(i,j)*x(j));
    #endfor
    #x(i) = (1/U(i,i))*(c(i) - otherTerms);
    x(i) = (1/U(i,i))*(c(i) - dot(U(i,i+1:r),x(i+1:r)));
  endfor
endfunction 

## Helper Function 
function x = solveEquation(A,b,method) 
  if (strcmp(method,"LR"))
    disp("Using LR Factorization with Backward and Forward Substitution");
    [L,U,P] = factorizationLR(A);
    c = forwardSubstitution(L,b(P));
    x = transpose(backwardSubstitution(U,c));
  elseif (strcmp(method,"naive"))
    disp("Using Guassian Elimination with Pivots");
    x = transpose(guassianElimination(A,b));
  endif
  
endfunction

## Helper Function 
function compareMethods(A,b) 
  start_lr = time();
  x_lr = solveEquation(A,b,"LR");
  end_lr = time();
  start_naive = time();
  x_naive = solveEquation(A,b,"naive");
  end_naive = time();
  
  duration_lr     = end_lr    - start_lr; 
  duration_naive  = end_naive - start_naive; 
  
  disp("Naive Method: ")
  disp(duration_naive)
  disp("Approximate Error: ")
  disp(sum(abs(b - (A*x_naive))));
  
  disp("")
  
  disp("LR Method: ")
  disp(duration_lr)
  disp("Approximate Error: ")
  disp(sum(abs(b - (A*x_lr))));
  
endfunction 
