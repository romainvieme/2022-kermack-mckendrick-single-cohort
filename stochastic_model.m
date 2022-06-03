% Please cite the following preprint (available at https://doi.org/10.48550/arXiv.2205.15634) if you use this code:
% A Kermack-McKendrick model with age of infection starting from a single or multiple cohorts of infected patients  
% by Jacques Demongeot, Quentin Griette, Yvon Maday, Pierre Magal
% arXiv e-print number 2205.15634. 2022
function [Time, IFlux] = stochastic_model(T, I, NS, tau, a0, beta0,beta01,AI,dt)
	t=0;

	% time step in days 
	l=size(I);
	m=1;
	m1=1;
	a1=0:dt:50; 
	a2=0:dt:50+dt; 
	beta1=beta(a1, a0, beta0,beta01);
	increment= dt;
	increment2 =0; 
	IFlux1=0; 
	while t<=T+1  

		% Next recovery event  

		k=find(I(2,:)>=0);
		I1=[I(1,k); I(2,k)];
		l=size(I1);          % we count the number of infected 

		Infected(2,m)=l(2); % we report the number of infected to plot it 
		Infected(1,m)=t;    % we report the instant at which this was observed 

		m=m+1;
		NS1=NS;
		increment=  increment-dt;
		p=1;

		if (l(2)>0)
			next_contact_rate = l(2)*NS*tau; % Old
			while l(2)>=1 && increment <=dt 

				% We count the average number of infectious




				next_contact = exprnd(1/next_contact_rate);

				i = randi([1 l(2)]);



				if rand()<= beta(I1(1,i), a0, beta0,beta01)  % Old
					% create a new infected invididual with age of infection 0
					I1(1,end+1) =0;
					I1(2,end+1)= exprnd(AI);

					% Update the Susceptible 
					NS1 = NS1-1;
					p=p+1; 
				end % else do nothing

				l=size(I1);


				%  next_contact_rate = Number_of_Infectious*NS*tau;
				next_contact_rate = l(2)*NS*tau;
				next_contact = exprnd(1/next_contact_rate);
				increment=increment+next_contact;
			end
		end
		if (increment2 <1)

			%IFlux(1,m1)=t; 

			IFlux1=(NS-NS1)+IFlux1; %the number of new infected in the interval of time [t,t+dt]
			increment2=increment2+dt; 
		else 
			Time(1,m1)=t; 
			IFlux(1,m1)=IFlux1;
			increment2=0;
			IFlux1=0; 
			m1=m1+1;
		end


		t=t+dt;
		I1(1,:)=I1(1,:)+dt;
		I1(2,:)=I1(2,:)-dt;
		I=I1; 
		NS=NS1;

	end
end
