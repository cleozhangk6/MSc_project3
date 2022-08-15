const express = require('express');  // require express
const { execFile } = require('child_process');  // child process module allows to execute any shell commands or scripts within nodejs
const upload = require('express-fileupload');  // module allows file upload
// const { response } = require('express');
const fs = require('fs');  // read file
const bodyParser = require('body-parser');
const download = require('download');

const app = express();  // create the app
var port = 3000
var host = '127.0.0.1'
app.listen(port, host, () => console.log(`Listening at port ${port}... Go to http://${host}:${port}/`));  //listen at a port (port number 3000)
app.use(express.static('public'));
app.use(upload())
app.use(express.json());
app.use(express.urlencoded({     // to support URL-encoded bodies
    extended: true
  }));

app.set('view engine', 'ejs');

app.get('/index.html', (req, res)=>{
    res.sendFile(__dirname+"/"+"index.html");
})


app.get('/search', (req, res) => {
    var uniprot = req.query.uniprot.toUpperCase();
    var filename = "AF-" + uniprot + "-F1-model_v2.pdb"
    // Url of the pdb file
    var url = 'https://alphafold.ebi.ac.uk/files/' + filename;
    // Path at which pdb will get downloaded
    var filePath = './public/uploads';
    download(url,filePath)
    .then(() => {
        console.log('Download from AlphaFold.org Completed');

        const parse = execFile('./run_parse_local.sh', [filename]);
                
        parse.stdout.on('data', (data) => {
            console.log(`stdout: ${data}`);
        });

        parse.stderr.on('data', (data) => {
            console.error(`stderr: ${data}`);
        });

        parse.on('close', (code) => {
            console.log(`child process exited with code ${code}`);                    
            fs.readFile(`public/outputs/domains_info.json`, 'utf8', (err, dataTxt) => {
                if (err) {
                    console.error(err);
                    return;
                }
                res.render('result', {data : {userQuery: uniprot,
                                        searchResults : JSON.parse(dataTxt)
                                        }})
            });
        }); 

})
})


app.post('/', (req, res) => {
    if (req.files) {
        console.log(req.files)
        var file = req.files.file
        var filename = file.name
        if (filename.startsWith("AF-")) {
            var uniprot = filename.split('-')[1]
        } else {
            var uniprot = "Unknown"
        }

        console.log(filename)

        file.mv('./public/uploads/'+filename, (err) => {
            if (err) {
                res.send(err)
            } else {
                const parse = execFile('./run_parse.sh', [filename]);
                
                parse.stdout.on('data', (data) => {
                    console.log(`stdout: ${data}`);
                });

                parse.stderr.on('data', (data) => {
                    console.error(`stderr: ${data}`);
                });
        
                parse.on('close', (code) => {
                    console.log(`child process exited with code ${code}`);                    
                    fs.readFile(`public/outputs/domains_info.json`, 'utf8', (err, dataTxt) => {
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

// // app.get('/:uniprot', (req, res) => {
// //     var uniprot = req.body.uniprot
// //     console.log(`hello ${uniprot}`)
// //     fs.readFile(`summary/${uniprot}_summary.txt`, 'utf8', (err, dataTxt) => {
// //         if (err) {
// //             console.error(err);
// //             return;
// //         }
// //         res.render('result', {data : {userQuery: uniprot,
// //                                 searchResults : dataTxt.split(/\r?\n/).slice(0,-1)}})
// //     });
// // })

// app.post('/', (req, res) => {
//     if (req.files) {
//         console.log(req.files)
//         var file = req.files.file
//         var filename = file.name
//         var uniprot = filename.split('-')[1]
//         console.log(filename)

//         file.mv('./uploads/'+filename, function (err) {
//             if (err) {
//                 res.send(err)
//             } else {
//                 // res.send("File Uploaded");
//                 // res.redirect("results.html");
//                 const parse = spawn('python3.9', ['temp/parse_mmcif_to_domains.py', `uploads/${filename}`]);

//                 parse.stdout.on('data', (data) => {
//                 console.log(`stdout: ${data}`);
//                 fs.readFile(`summary/${uniprot}_summary.txt`, 'utf8', (err, dataTxt) => {
//                     if (err) {
//                         console.error(err);
//                         return;
//                     }
//                     // console.log(dataTxt);

//                     res.render('result', {data : {userQuery: uniprot,
//                                             searchResults : dataTxt.split(/\r?\n/).slice(0,-1)}})
//                 });
//                 });
        
//                 parse.stderr.on('data', (data) => {
//                 console.error(`stderr: ${data}`);
//                 });
        
//                 parse.on('close', (code) => {
//                 console.log(`child process exited with code ${code}`);
//                 }); 
                
//             }
//         })


//     }
// })
