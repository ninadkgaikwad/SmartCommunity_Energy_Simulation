function [] = SendEmail_From_MATLAB(EmailInput_Struct)

% Author: Ninad Kiran Gaikwad
% Date: May/21/2020
% Description: SendEmail_From_MATLAB - Send Email from MATLAB

%% Sending Email from MATLAB

%% Getting Required Variables from Input Struct

% Email Sender Information
Email_Sender = EmailInput_Struct.Email_Sender;
Email_Password_Sender = EmailInput_Struct.Email_Password_Sender;
SMTP_Server_Sender=EmailInput_Struct.SMTP_Server_Sender;

% Email Contents
Email_ReceiverList_Cell=EmailInput_Struct.Email_ReceiverList_Cell;
Email_Subject=EmailInput_Struct.Email_Subject;
Email_Text_Vector=EmailInput_Struct.Email_Text_Vector;
Email_Attachment_Cell=EmailInput_Struct.Email_Attachment_Cell;

%% Sending Email from MATLAB

% Set Internet Properties
setpref('Internet','SMTP_Server',SMTP_Server_Sender);
setpref('Internet','E_mail',Email_Sender);
setpref('Internet','SMTP_Username',Email_Sender);
setpref('Internet','SMTP_Password',Email_Password_Sender);

% Java Properties to be set
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send the email.  Note that the first input is the address you are sending the email to
sendmail(Email_ReceiverList_Cell,Email_Subject,Email_Text_Vector,Email_Attachment_Cell)
end

