const express = require('express');  // require express
const { execFile } = require('child_process');  // child process module allows to execute any shell commands or scripts within nodejs
const upload = require('express-fileupload');  // module allows file upload
const { response } = require('express');
const fs = require('fs');  // read file

const app = express();  // create the app
var port = 20236
var host = '0.0.0.0'
app.listen(port, host, () => console.log(`Listening at port ${port}... Go to http://${host}:${port}/`));  //listen at a port (port number 3000)
app.use(express.static('public'));
app.use(upload())
app.use(express.json())

// app.get('/', (req, res) => {
//     res.sendFile(__dirname + 'index.html')
// })

app.set('view engine', 'ejs');
// app.get('/:userQuery', (req, res)=>{
//     fs.readFile(`summary/${req.params.userQuery}_summary.txt`, 'utf8', (err, dataTxt) => {
//         if (err) {
//             console.error(err);
//             return;
//         }
//         // console.log(dataTxt);
//         res.render('result', {data : {userQuery: req.params.userQuery,
//                                 searchResults : dataTxt.split(/\r?\n/).slice(0,-1)}})
//     });
// })


app.post('/', async (req, res) => {
    if (req.files) {
        console.log(req.files)
        var file = req.files.file
        var filename = await file.name
        var uniprot = await filename.split('-')[1]
        console.log(filename)

        await file.mv('./uploads/'+filename,  async (err) => {
            if (err) {
                res.send(err)
            } else {
                // res.send("File Uploaded");
                // res.redirect("results.html");
                const parse = execFile('./run_parse.sh', [filename]);
                
                parse.stdout.on('data', async (data) => {
                    console.log(`stdout: ${data}`);
                });

                parse.stderr.on('data', (data) => {
                    console.error(`stderr: ${data}`);
                });
        
                parse.on('close', async (code) => {
                    console.log(`child process exited with code ${code}`);                    
                    fs.readFile(`outputs/domains_info.json`, 'utf8', (err, dataTxt) => {
                        if (err) {
                            console.error(err);
                            return;
                        }
                        // const dataJson = JSON.parse(dataTxt)
                        res.render('result', {data : {userQuery: uniprot,
                                                searchResults : JSON.parse(dataTxt)
                                                }})
                    });
                }); 
                

                
            }
        })


    }
})

