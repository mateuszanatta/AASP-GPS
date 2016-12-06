function data = orderStruct(data)
% Order satellite and channel structs in PRN order, making comparison
% easier
%
%[] = orderStruct(data)

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Written by Alexandre Serio Buscher,   IME, Rio de Janeiro Brazil
%                                       TU Ilmenau, Germany
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

tam = length(data);
% selection sort
for i = 1:(tam-1)
   min = data(i).PRN;
   for j = (i+1):tam
       if (data(j).PRN < min  &&  data(j).PRN ~= 0 )
           min = data(j).PRN;
           index = j;
       end %if (data(j).PRN > min)
   end %for j = (i+1):tam
   if (data(i).PRN ~= min  )
       temp = data(i);
       data(i) = data(index);
       data(index) = temp;
   elseif data(i).PRN == 0 
       break;
   end %if (data(i).PRN ~= min)
end %for i = 1:tam 
