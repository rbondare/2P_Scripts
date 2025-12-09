%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   LCr class for communicating LightCrafter (Texas Instrument) via TCP/IP
%
%   ----- Usage example -----
%   >> l = LCr;
%   >> l.connect; % make sure to "disconnect" LCr on GUI by TI before connecting via matlab
%   >> l.setDisplayMode('02'); % HDMI external video mode
%
%   >> l.getImageProcessingStatus; % check the image processing status
%   >> l.disableImageProcessing;
%
%   >> l.getImageProcessingStatus; % make sure that the image processing is all disabled!!
%
%   >> l.getLEDCurrent; % LED current status
%   >> l.setLEDCurrent(1, 274, 274); % set LED current (0-274)
%   
%   >> l.disconnect; % make sure to disconnect LCr before using GUI provided by TI
%
%   ----- references -----
%       DLP LightCrafter DM365 Command Interface Guide
%       DLPC300 Programmers Guide
%       http://www.ti.com/
%
%   ----- copy right -----
%   Original author: Jan Winter, TU Berlin, FG Lichttechnik, j.winter@tu-berlin.de
%   Modified by: Hiroki Asari, Markus Meister Lab, Caltech, hasai@caltech.edu

classdef LCr_matlab < handle
    
    properties %(Hidden)
        tcpConnection
    end
    
    methods

        %constructor
        function obj = LCr_matlab()
        end

        function connect( obj )
%             echotcpip('on', 21845); % start tcpip server
            obj.tcpConnection = tcpip( '192.168.1.100', 21845 );
%             obj.tcpConnection.BytesAvailableFcn = @instrcallback
%             obj.tcpConnection.BytesAvailableFcnCount = 7;
%             obj.tcpConnection.BytesAvailableFcnMode = 'byte';
            fopen( obj.tcpConnection );
        end

        function disconnect( obj )
            fclose( obj.tcpConnection );
%             echotcpip('off'); % stop tcpip server
        end

        function header = createHeader( obj )
            header = uint8( zeros( 6, 1 ) );
        end

        function modifiedHeader = computePayloadLength( obj, header, payload )
            header( 5 ) = uint8( mod( length( payload ), 256 ) ); %payloadLength LSB
            header( 6 ) = uint8( floor( length( payload ) / 256 ) ); %payloadLength MSB
            modifiedHeader = header;
        end

        function modifiedPacket = appendChecksum( obj, packet )
            checksum = mod( sum( packet ), 256 );
            modifiedPacket = [ packet; checksum ];
        end

        function getVersion( obj, Version )

            if (~ischar(Version) && (length(Version)~=2))
                disp('Version must be a 2 digit hex string 00, 10, or 20');
                return;
            end

            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '04' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '00' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( hex2dec( Version ) ); %payload
            header = obj.computePayloadLength( header, payload);
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            
            obj.getData( packet );
        end

        function setDisplayModeStatic( obj )
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '01' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
            header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
            header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( hex2dec( '00' ) ); %payload
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            obj.sendData( packet );

        end

        function setDisplayModeInternalPattern( obj )
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '01' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
            header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
            header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( hex2dec( '01' ) ); %payload
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            obj.sendData( packet );

        end
        
        function setDisplayMode( obj, DisplayMode )

            if (~ischar(DisplayMode) && (length(DisplayMode)~=2))
                disp('Display mode must be a 2 digit hex string in range 00 to 04');
                return;
            end

            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '01' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( hex2dec( DisplayMode ) ); %payload
            header = obj.computePayloadLength( header, payload);
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            obj.sendData( packet );
        end
        
        function getDisplayMode( obj )
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '04' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '01' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '00' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            packet = obj.appendChecksum( header ); %packet
            obj.getData( packet );
        end

        function setLEDCurrent( obj, R, G, B )
            
            % 0x0 to 0x112 (w/o active cooling) to 0x400 (w/ active cooling) 
            if ischar(R), R=hex2dec(R);end 
            if R<0, R=0; disp('minimum Red LED curent is 0.');end
            if R>274 % R=274; 
                warning('Red LED current %u is above max recommended range of 274 without active cooling',R);end
            R = dec2hex(R,3);
            
            if ischar(G), G=hex2dec(G);end 
            if G<0, G=0; disp('minimum Green LED curent is 0.');end
            if G>274 % G=274; 
                 warning('Green LED current %u is above max recommended range of 274 without active cooling',G);end
            G = dec2hex(G,3);
            
            if ischar(B), B=hex2dec(B);end 
            if B<0, B=0; disp('minimum Blue LED curent is 0.');end
            if B>274 % B=274; 
                 warning('Blue LED current %u is above max recommended range of 274 without active cooling',B);end
            B = dec2hex(B,3);
            
%             if (~ischar(redCurrent) && (length(redCurrent)~=3))
%                 disp('Red LED curent must be 000 to 400 (1024 dec)');
%                 return;
%             end
%             if (~ischar(greenCurrent) && (length(greenCurrent)~=3))
%                 disp('Green LED curent must be 000 to 400 (1024 dec)');
%                 return;
%             end
%             if (~ischar(blueCurrent) && (length(blueCurrent)~=3))
%                 disp('Green LED curent must be 000 to 400 (1024 dec)');
%                 return;
%             end

            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '04' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '06' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( [ hex2dec( R(2:3) ); hex2dec( R(1) ); hex2dec( G(2:3) ); hex2dec( G(1) ); hex2dec( B(2:3) ); hex2dec( B(1) ) ] ); %payload
            header = obj.computePayloadLength( header, payload);
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            obj.sendData( packet );

        end

        function getLEDCurrent( obj )
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '04' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '04' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '00' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            packet = obj.appendChecksum( header );
            obj.getData( packet );
        end        
        
        function setInternalPattern( obj, pattern )

%             obj.setDisplayModeInternalPattern();
            obj.setDisplayMode( '01' );
            
            if ( ~ischar( pattern ) && ( length(pattern) ~= 2 ) )
                disp('pattern must be a 2 digit hex string in range 00 to 0D')
                return;
            end

            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '03' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
%             header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
%             header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( hex2dec( pattern ) ); %payload
            header = obj.computePayloadLength( header, payload);
            packet = obj.appendChecksum( [ header; payload ] ); %packet
            obj.sendData( packet );

        end

        function setStaticColor( obj, R, G, B )

%             obj.setDisplayModeStatic();
            obj.setDisplayMode( '00' );
            
            if ( ~ischar( R ) && ( length(R) ~= 2 ) )
                disp('R must be a 2 digit hex string in range 00 to FF');
                return;
            end
            if ( ~ischar( G ) && ( length(G) ~= 2 ) )
                disp('G must be a 2 digit hex string in range 00 to FF');
                return;
            end
            if ( ~ischar( B ) && ( length(B) ~= 2 ) )
                disp('B must be a 2 digit hex string in range 00 to FF');
                return;
            end

            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '06' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
            header( 5 ) = uint8( hex2dec( '04' ) ); %payloadLength LSB
            header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            payload = uint8( [ hex2dec( B ); hex2dec( G ); hex2dec( R ); hex2dec( '00' ) ] ); %payload         
            packet = obj.appendChecksum( [ header; payload ] );
            %packet
            obj.sendData( packet );

        end

        function setBMPImage( obj, imageData )

%             obj.setDisplayModeStatic( );
            obj.setDisplayMode( '00' );

            MAX_PAYLOAD_SIZE = 65535;
            numberOfChunks = ceil( length( imageData ) / 65535 );
            chunkArray = cell( numberOfChunks, 1 );
            for i = 1 : numberOfChunks
                currentLength = length( imageData );
                if( currentLength > MAX_PAYLOAD_SIZE )
                    chunkArray{ i } = imageData( 1 : MAX_PAYLOAD_SIZE );
                    imageData = imageData( MAX_PAYLOAD_SIZE + 1 : end );
                else
                    chunkArray{ i } = imageData( 1 : end );
                end
            end

            for currentChunkIndex = 1 : numberOfChunks

                currentChunk = chunkArray{ currentChunkIndex };

                header = obj.createHeader();
                header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
                header( 2 ) = uint8( hex2dec( '01' ) ); %CMD1
                header( 3 ) = uint8( hex2dec( '05' ) ); %CMD2
                header = obj.computePayloadLength( header, currentChunk );

                %append flag
                if( numberOfChunks == 1 )
                    header( 4 ) = uint8( hex2dec( '00' ) ); %flags
                else
                    if( currentChunkIndex == 1 )
                        disp('FIRST CHUNK')
                        header( 4 ) = uint8( hex2dec( '01' ) ); %flags
                    elseif( currentChunkIndex == numberOfChunks )
                        disp('LAST CHUNK')
                        header( 4 ) = uint8( hex2dec( '03' ) ); %flags
                    else
                        disp('OTHER CHUNK')
                        header( 4 ) = uint8( hex2dec( '02' ) ); %flags
                    end
                end

                packet = obj.appendChecksum( [ header; currentChunk ] );
                obj.sendData( packet );
            end

        end

        function getImageProcessingStatus( obj )
            
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '04' ) );	%packet type
            header( 2 ) = uint8( hex2dec( 'FF' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '00' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
            header( 5 ) = uint8( hex2dec( '01' ) ); %payloadLength LSB
            header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            
            packet = '50'; % automatic gain control
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            fwrite( obj.tcpConnection, packet ) ; pause(1)
            if obj.tcpConnection.BytesAvailable>0,
                LCr_out = fread(obj.tcpConnection, obj.tcpConnection.BytesAvailable);
                if LCr_out(1)==5 && LCr_out(7)==6, 
                    disp('automatic gain control : disabled');
                else 
                    disp('automatic gain control : enabled');
                end
            end
            
            packet = '7E'; % temporal dithering
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            fwrite( obj.tcpConnection, packet ) ; pause(1)
            if obj.tcpConnection.BytesAvailable>0,
                LCr_out = fread(obj.tcpConnection, obj.tcpConnection.BytesAvailable);
                if LCr_out(1)==5 && LCr_out(7)==2, 
                    disp('temporal dithering : disabled');
                else 
                    disp('temporal dithering : enabled');
                end
            end
            
            packet = '5E'; % color coordinate adjustment
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            fwrite( obj.tcpConnection, packet ) ; pause(1)
            if obj.tcpConnection.BytesAvailable>0,
                LCr_out = fread(obj.tcpConnection, obj.tcpConnection.BytesAvailable);
                if LCr_out(1)==5 && LCr_out(7)==0, 
                    disp('color coordinate adjustment : disabled');
                else 
                    disp('color coordinate adjustment : enabled');
                end
            end
        end   
        
        function disableImageProcessing( obj )
            
            header = obj.createHeader();
            header( 1 ) = uint8( hex2dec( '02' ) );	%packet type
            header( 2 ) = uint8( hex2dec( 'FF' ) ); %CMD1
            header( 3 ) = uint8( hex2dec( '00' ) ); %CMD2
            header( 4 ) = uint8( hex2dec( '00' ) ); %flags
            header( 5 ) = uint8( hex2dec( '05' ) ); %payloadLength LSB
            header( 6 ) = uint8( hex2dec( '00' ) ); %payloadLength MSB
            
            packet = ['50'; '06'; '00'; '00'; '00']; % automatic gain control
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            obj.sendData( packet );
            
            packet = ['7E'; '02'; '00'; '00'; '00']; % temporal dithering
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            obj.sendData( packet );
            
            packet = ['5E'; '00'; '00'; '00'; '00']; % color coordinate adjustment
            packet = obj.appendChecksum( [ header; uint8( hex2dec( packet ) ) ] );
            obj.sendData( packet );
        end  
        
        function sendData( obj, packet )

            MAX_SIZE = 512; %limit packet size
            buffer = packet;
            while (~isnan(buffer))
                if( length(buffer) > MAX_SIZE )
                    currentPacket = buffer( 1 : MAX_SIZE );
                    buffer = buffer( MAX_SIZE + 1 : end );
                else
                    currentPacket = buffer( 1 : end );
                    buffer = NaN;
                end
                fwrite( obj.tcpConnection, currentPacket ); pause(1);
                if obj.tcpConnection.BytesAvailable>0,
                    LCr_out = fread(obj.tcpConnection, obj.tcpConnection.BytesAvailable);
                    if LCr_out(1)==3, 
%                         disp('data written onto LCr:');
%                         disp( currentPacket );
                    end
                end
            end
        end
        
        function getData( obj, packet )
            
            MAX_SIZE = 512; %limit packet size
            buffer = packet;
            while (~isnan(buffer))
                if( length(buffer) > MAX_SIZE )
                    currentPacket = buffer( 1 : MAX_SIZE );
                    buffer = buffer( MAX_SIZE + 1 : end );
                else
                    currentPacket = buffer( 1 : end );
                    buffer = NaN;
                end
                fwrite( obj.tcpConnection, currentPacket ); pause(1);
                if obj.tcpConnection.BytesAvailable>0,
                    LCr_out = fread(obj.tcpConnection, obj.tcpConnection.BytesAvailable);
                    if LCr_out(1)==5,
                        disp('data read from LCr (header):');
                        disp(LCr_out(1:6));
                        disp('data read from LCr (data):');
                        disp(LCr_out(7:end-1));
                        disp('data read from LCr (checksum):');
                        disp(LCr_out(end));
                    end
                end
            end
        end
        
    end % methods

end % classdef