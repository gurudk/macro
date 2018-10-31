%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is the Matlab code for solving the stochastic growth model     %
%  with log utility, Cobb-Douglas production function,                 %
%  incomplete depreciation, and a Markov chain process                 %
%  for productivity shocks.                                            %
%  The code implements the standard value function iteration method.   %
%                                                                      %             
%   Author  : Zhang Yongji                                             %
%   History : 2017-06-25  modified from Kevin X.D. Huang's loops1.m    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%-----------------------Main ideas and steps-------------------------%%
%  The bellman equation:                                               %
%           V(s,k) = max(u(Ct)) + ¦ÂEs'|s(V(s',k'))                     %
%                     k'                                               %
%                                                                      %
%  Main idea for value function iterator;                              %
%   First, choose an endogenous state space for capital, begin the     %
%   iterator for V(k)=0,and find a sequence of k',and compuate V1(k)   %
%   for next iterator.Over and over, stop util that the distance of    %
%   lastest two iterators is smaller than the error we specified before.
%                                                                      %
%  Main steps for this program:                                        %
%   Step 1: Set up for init parameters.                                %
%   Step 2: Value function iterate                                     %
%   Step 3: Draw the results of value function and policy function     %
%%--------------------------------------------------------------------%%

% It's a good habit to clear variable value before run a matlab script
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%--Step 1--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Global constants
epsilon_    = -10^10;               % very very small value

%Set parameter values
alpha       = 0.3;                  % production parameter
beta        = 0.98;                 % subjective discount factor 
delta       = 0.025;                % capital depreciation rate

%Exogenous state space and markov transition matrix
s           = [-3 3];
p(1:2,1:2)  = 0.5;

%Steady state capital
kbar        =(1/(alpha*beta)-(1-delta)/alpha)^(1/(alpha-1));

%Choose endogenous state space
mink        = 0.1*kbar;
maxk        = 10*kbar;   
nk          = 5;                        % generate nk x nk grid
gridk       = linspace(mink,maxk,nk);   % grid generation

%Set iteration criterion
maxit       = 10^5;                      %maximum for iterator
tol         = 10^-5;                     %tolerence for iterator

%Initial value function
value_function(2,nk,1)                      = 0;
policy_function(2,nk,1)                     = 0;
lastest_two_iterator_results(1:2,1:2*nk)    = 0;

%Value function iteration
it          = 1;
diff        = 10^10;


%%%%%%%%%%%%%%%%%%%%%%%%%--Step 2--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while it<=maxit & diff>tol;
      it = it+1;
      %For every state in state space
      for j = 1 : length(s);
          %Traverse the grid for state j
          for k_current=1:nk;
              %Init value function at begin.
              value_function(j,k_current,it) = epsilon_;
              
              %Find max value of all k_next, and record the result of
              %value function and policy function
              for k_next=1:nk;
                  %Compute consume for current period
                  cons = exp(s(j))*gridk(k_current)^alpha+(1-delta)*gridk(k_current)-gridk(k_next);
                  if (cons > 0)
                    %Compute Vk given k_next
                    vk = log(cons)+beta*(value_function(:,k_next,it-1)'*p(:,j));
                    if vk >= value_function(j,k_current,it)
                        %Find the max Vk for k_current!!!
                        %Just record value function and policy funcion
                        value_function(j,k_current,it) = vk;
                        policy_function(j,k_current,it) = gridk(k_next);
                    end;
                  end; % end of cons if
              end;% end of k_next
          end;% end of k_current
      end;% end of state s
      
      %Compute indicator for this round
      indicator = mod(it,2);
      if (indicator == 0)
        indicator = 2;
      end
      
      %Vectorization of this value function result for this round
      lastest_two_iterator_results(indicator,:) = reshape(value_function(:,:,it),[1,2*nk]);
      
      %Compute difference of lastest two iterator results
      diff = max(abs(lastest_two_iterator_results(1,:)-lastest_two_iterator_results(2,:)));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%--Step 3--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot figure
k_current=1:1:nk;
subplot(2,1,1);
plot(gridk(k_current),value_function(1,k_current,4),'-.',gridk(k_current),value_function(1,k_current,5),'-',gridk(k_current),value_function(2,k_current,4),'-o',gridk(k_current),value_function(2,k_current,5),'-*');
xlabel('Capital Stock k Today')
ylabel('Value Function')
title('Value Function')

subplot(2,1,2);
plot(gridk(k_current),policy_function(1,k_current,4),'-.',gridk(k_current),policy_function(1,k_current,5),'-',gridk(k_current),policy_function(2,k_current,4),'-o',gridk(k_current),policy_function(2,k_current,5),'-*');
xlabel('Capital Stock k Today')
ylabel('Policy Function')
title('Policy Function')
