% SynScan Pro TCP/IP 설정
host = '192.168.4.1'; % SynScan Pro의 IP 주소
port = 11881;         % SynScan Pro의 포트 번호 (기본값)

% TCP/IP 연결 설정
client = tcpclient(host, port);


% RA 설정 (예: 10:00:00)
raCommand = ':Sr10:00:00#';
sendCommand(client, raCommand);

% DEC 설정 (예: +20:00:00)
decCommand = ':Sd+20:00:00#';
sendCommand(client, decCommand);

% GOTO 명령 실행
gotoCommand = ':MS#';
sendCommand(client, gotoCommand);

% 연결 닫기
clear client;
disp('GOTO 명령 완료');



% 명령 보내는 함수 정의
function sendCommand(client, command)
    write(client, uint8(command), 'uint8'); % 명령 전송
    pause(0.1); % SynScan Pro가 명령을 처리할 시간 대기
    response = read(client); % 응답 읽기
    disp(['Response: ', char(response')]); % 응답 출력
end